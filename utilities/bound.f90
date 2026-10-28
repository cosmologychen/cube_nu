module bound_test
   ! 生产求解复制自 halo_bound_check_handoff/review/bound_test_diagnostic_review.f90。
   ! 仅移除诊断、串行参照、旧链表兼容及通用测试输入验证，数值计算语句不变。
   ! 调用方提供 xv、pid、成员CRS及快照 a_bound/b_bound；坐标须已展开周期边界。
   ! 原始粒子相空间、成员CRS和稳定ID在求解中不改写。
   ! 物理 CPU affinity 由调用方运行配置负责；不在本模块内绑定 CPU。
   use parameters, only: ng,ng_global,nn,box,ratio_cs,h0,np_halo_min,pi,app,n_int
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use omp_lib, only: omp_get_thread_num, omp_get_num_threads, omp_get_wtime, &
      omp_set_max_active_levels, omp_in_parallel
   use iso_c_binding, only: c_loc, c_f_pointer
   implicit none

   ! worker 是一个串行处理单 Halo 的 OpenMP 线程；16 不是每次 FFT 的线程数。
   integer,parameter :: fofcore=64
   ! 数组顺序为 [48,64,128,256]；默认 [6,3,5,2] 对应文档的 256/128/64/48=2/5/3/6。
   integer :: bound_workers(4)=[6,3,5,2]
   ! s 为全队线程编号+1；worker_group(s)=网格组，worker_slot(s)=该组内独占槽号。
   integer :: worker_group(fofcore)=0,worker_slot(fofcore)=0
   ! 原位 FFT：rho 的第一维 n+2 为实数 padding；ft 只是同一内存的复数视图，不另占一份。
   type bound_fft_storage
      real(4),pointer,contiguous :: rho(:,:,:,:)=>null()
      complex(4),pointer,contiguous :: ft(:,:,:,:)=>null()
   endtype
   type(bound_fft_storage) :: bound_storage(4)
   ! [阶段 2.1 生产接入 / 阶段 2.2 动态去束缚]: P3M 核分裂截断距离常数 h_split = 3.5 格点
   real(8),parameter :: bound_h_split=3.5d0
   private :: bound_c2_delta_depth, &
      bound_fft_wave_number, bound_assignment_window_sq

   !----------------------------------------------------------------------
   ! 动态多分辨率孤立 PM 网格方案（替代旧的 ng/[16,8,6,4] 固定分档）
   !
   ! 1. 候选孤立 FFT 工作网格每维节点数固定为 48^3 / 64^3 / 128^3 / 256^3；
   !    只启用满足 n_grid <= ng 的网格，最大网格 n_max = min(256, ng)。
   !    n_grid>ng 的类别由 ngrid 禁用；批入口另外要求 ng>=256，
   !    不为禁用类别创建 FFT 计划或分配 PM 工作区。
   !    bn 的 min() 仅限制尺寸元数据，不表示批入口支持所有较小 ng。
   ! 2. D_max = 15 Mpc 只限制 PM 候选的最大轴向跨度；直接求和不受此覆盖限制。
   !    15 Mpc 对应的标准格点跨度 S_15 = 15 / L_ng = 15 * h0 * ng_global / box。
   !    此处 box 以 Mpc/h 给出，L_ng 与 15 Mpc 均取共动 Mpc，两端约定一致。
   ! 3. 最大网格有效间隔数 m_max = n_max/2 - 3，最大（终止）间距
   !    spacing_max = S_15 / m_max，使最大档有效覆盖严格等于 15 Mpc。
   ! 4. 各网格档位由 bound_grid_levels 按几何序列 1/32, 1/16, 1/8, 1/4, 1/2, ...
   !    生成：只保留物理覆盖严格小于 15 Mpc 的档，最大启用网格再追加精确
   !    spacing_max。所有档位按物理覆盖升序组成全局阶梯，相邻档连续覆盖。
   !
   ! 生产外层网格上限为256，不实现或分配512外层网格。
   !----------------------------------------------------------------------
   integer,parameter :: cand(4) = [48, 64, 128, 256]
   integer,parameter :: ngi = int(ng)
   ! 实际启用的网格档数：ng>=256 为 4，256>ng>=128 为 3（要求 ng>=48）
   integer,parameter :: ngrid = 1 + merge(1,0,ngi>=64) + merge(1,0,ngi>=128) + merge(1,0,ngi>=256)
   ! 各类别尺寸元数据：min(候选尺寸, ng)；可写工作区在批入口按实际组槽数分配。
   integer,parameter :: bn(4) = min(cand, ngi)
   ! 最大启用网格每维节点数、最大有效间隔数、15 Mpc 标准格点跨度及终止间距
   integer,parameter :: n_max = min(256, ngi)
   integer,parameter :: m_max = n_max/2 - 3
   real(8),parameter :: s_15 = 15d0 * dble(h0) * dble(ng_global) / dble(box)
   real(8),parameter :: spacing_max = s_15 / dble(m_max)

   ! 每类仅分配该组槽数；复数视图与实数 padding 工作区共享原位 FFT 内存。

   !----------------------------------------------------------------------
   ! 数据结构 1: bound_kernel_type（网格引力核与自身势基核）
   ! 存放某一网格类别预计算出的频域 Green 函数和自身势基核。
   ! 该结构在初始化后完全只读，同类别的所有工作线程/槽位安全共享，无需重复分配。
   !----------------------------------------------------------------------
   type bound_kernel_type
      integer :: n=0             ! 孤立 PM 网格每维总节点数（含后半域零填充空间，如 48, 64 等偶数）
      integer :: nlevel=0        ! 该网格类别按动态网格规则启用的 spacing 档位数（随 ng 与 box 变化）
      integer :: ng_standard=0   ! 核初始化时的单节点 ng；用于现有缓存复用检查
      real(8),allocatable :: spacing(:)       ! 各 level 的网格间距比例（以标准格点间距为基准）
      real(8) :: softening=0d0   ! 物理软化长度（标准格点单位 app）
      ! C1 势最优频域 Green 函数核与实空间 5x5x5 自项基核 (TSC order=2, h_split=3.5)
      complex(4),allocatable :: green_phi(:,:,:)
      real(4) :: self_phi(-2:2, -2:2, -2:2) = 0.0
   endtype

   !----------------------------------------------------------------------
   ! 数据结构 2: bound_mesh_type（工作槽绑定的网格元数据与 FFT 句柄）
   ! 每个并发工作槽位（或 OpenMP 线程）拥有一套 FFT 计划与状态描述。
   ! rho 与 ft 是同一原位 FFT 内存的实数/复数视图，工作槽按组内编号取连续切片。
   ! FFT 计划绑定各槽的可写工作区；Green 和 self_phi 则按类别只读共享。
   !----------------------------------------------------------------------
   type bound_mesh_type
      integer(8) :: plan=0,iplan=0              ! 绑定的单精度 FFTW 正向 (r2c) 与逆向 (c2r) 变换句柄计划
      integer :: grid_id=0                      ! 槽的 FFT 类别；任务方法另由 bound_mem%method_gid 表示
      integer :: level=0                        ! 当前类别合法 spacing 表中的档位序号
      integer :: required_grid=0                ! 旧接口保留字段；当前固定批入口不重新调度或迁移任务
      real(8) :: spacing=0d0                    ! 当前选定 level 的真实网格间距比例
      real(8) :: span=0d0                       ! 当前候选 Halo 在展开坐标下的三维包围盒最大边长（标准格点单位）
   endtype

   !----------------------------------------------------------------------
   ! 数据结构 3: bound_member_type（工作槽独占的粒子成员缓冲）
   ! 存放当前候选 Halo 的成员索引及束缚状态掩码，独立于网格尺寸。
   !----------------------------------------------------------------------
   type bound_member_type
      integer,allocatable :: member(:)    ! 候选 Halo 成员在原相空间全局数组 xv 中的原始粒子下标
      logical,allocatable :: keep(:)      ! 束缚状态掩码
      integer,allocatable :: order(:),first_tsc(:,:),just_del_list(:)
      real(4),allocatable :: q(:,:),self(:),phi(:),energy(:),wt_tsc(:,:,:),phi_c1(:),phi_c2(:)
      ! phi_direct 跨外层轮保留双精度势；能量判定用的 phi 是它的单精度副本。
      real(8),allocatable :: phi_direct(:)
      real(4),allocatable :: phi_c2_team(:,:)
      ! anchor 是原快照坐标参考点；origin 是相对 anchor 的标准格坐标下 PM 源域原点。
      real(4) :: anchor(3)=0.0,origin(3)=0.0
      integer :: method_gid=0
   endtype

   !----------------------------------------------------------------------
   ! [C2串行优化五] 数据结构 3b: bound_cell_type（工作槽独占的连续 cell/CSR 结构）
   ! 每 Halo 在 spacing 确定后构建一次, 去束缚全程粒子位置不变, 跨 pass 复用;
   ! 替代 hoc/ll 传统链表, 提供连续内存访问与缓存局部性。
   ! 成员表先按稳定 ID 排序；cell 存成员局部索引，索引本身不是全局稳定 ID。
   !----------------------------------------------------------------------
   type bound_cell_type
      integer :: ncell(3) = 0              ! 每维 cell 数
      integer :: ncell_total = 0           ! cell 总数
      real(8) :: u_min(3) = 0d0            ! cell 网格原点 (缩放坐标 u = q/spacing)
      integer,allocatable :: cell_id(:)        ! (np_fof) 每粒子的 1D cell 编号
      integer,allocatable :: cell_start(:)     ! (ncell_total+1) CSR 起始偏移，最后一项为结束位置加一
      integer,allocatable :: cell_particles(:) ! (np_fof) cell 内按已排序成员的局部索引升序
      integer :: ntasks=0
      integer,allocatable :: task_c1(:),task_p1(:),task_p2(:)
   endtype
   ! [C2审核修正] 工作槽数组改为编译期固定长度 (fofcore 为 parameter):
   ! 外层数组随模块静态创建, 消除并行冷启动路径上的共享 allocate 竞争;
   ! 每槽 allocatable 组件由槽位独占线程管理, 同一槽位不存在并发修改。
   type(bound_cell_type) :: bound_cells(fofcore)
   ! [C2串行优化五]: 13 半邻域 stencil (与链表版一致, 提为模块级供 cell 版函数共用)
   integer,parameter :: bound_c2_stencil(3, 13) = reshape([ &
      -1,-1, 1,   0,-1, 1,   1,-1, 1, &
      -1, 0, 1,   0, 0, 1,   1, 0, 1, &
      -1, 1, 1,   0, 1, 1,   1, 1, 1, &
      -1, 1, 0,   0, 1, 0,   1, 1, 0, &
      1, 0, 0 ], [3, 13])

   !----------------------------------------------------------------------
   ! 数据结构 4: bound_schedule_type（多 Halo 任务调度表）
   ! 记录一次性分类结果与固定组内队列；保留 slot/nwave 字段，但当前只有一个固定执行批。
   !----------------------------------------------------------------------
   type bound_schedule_type
      integer,allocatable :: group(:)      ! 每个候选匹配的调度工作组 (1..ngrid，0 为未分配/异常；直接求和归 1 组)
      integer,allocatable :: method_gid(:) ! 每个候选选定的求势方法/网格类别 (0=直接求和, 1..ngrid=PM 网格类别)
      integer,allocatable :: level(:)      ! 每个候选选定的档位序号 (直接求和为 0)
      real(8),allocatable :: spacing(:)    ! 实际单精度 spacing 的双精度镜像
      real(4),allocatable :: anchor(:,:),origin(:,:)
      ! capacity 含本组所有任务；pm_capacity 只含 PM 任务，直接任务无需 TSC/C1/C2 缓冲。
      integer :: pm_capacity(4)=0
      integer,allocatable :: flag(:)       ! 旧接口保留字段；当前批入口要求初始候选已满足门槛
      ! 跨度超过 15 Mpc 不再用标志跳过，而是打印并非零终止
      integer,allocatable :: queue(:)   ! 各组按估算工作量降序排列的候选索引队列
      integer,allocatable :: slot(:,:)  ! 当前仅 slot(:,1)：每个 worker 固定所属的网格工作组
      real(4),allocatable :: span(:)    ! 各候选实际展开后的最大包围盒跨度（标准格点单位）
      real(8),allocatable :: cost(:)    ! 各候选估算的工作量分数
      integer :: first(4)=1,last(4)=0   ! 各组在 queue 中的起止指针
      integer :: capacity(4)=0          ! 各组所需最大候选成员数（供 bound_init_rho 分配缓存）
      integer :: nworker(4)=0           ! 各组分配的并发线程数
      integer :: nwave=0                ! 旧接口保留字段；当前固定批设置为 1
      real(8) :: work(4)=0d0            ! 各组总工作量分数
   endtype

   ! 模块级全局工作区与物理上下文变量：
   ! bound_gk: 4 类网格各自的预计算引力核（只读共享）
   type(bound_kernel_type) :: bound_gk(4)
   ! bound_mesh: 网格工作槽数组，大小为 fofcore，各槽独占使用
   type(bound_mesh_type),allocatable :: bound_mesh(:)
   ! bound_mem: 成员粒子工作槽数组，大小为 fofcore，各槽独占使用
   type(bound_member_type),allocatable :: bound_mem(:)
   ! a_bound: 当前快照尺度因子 a；b_bound: 哈勃相对速度流修正系数 B = H(a) * a^2 * T0
   real(8) :: a_bound=1d0,b_bound=0d0

   !----------------------------------------------------------------------
   ! 全局档位阶梯（按物理覆盖升序，由 bound_build_ladder 惰性建立一次）
   ! 每一行是一档：网格类别、该类别内部的 level 序号、spacing 与物理覆盖
   ! （覆盖 = 该网格有效间隔数 × spacing，标准格点单位）。调度与档位选择
   ! 都读取同一张表，保证"选出的网格"与"槽位内的档位"严格一致。
   !----------------------------------------------------------------------
   integer,parameter :: nlev_max = 64
   integer :: bound_ntotal = 0
   integer :: bound_tab_gid(nlev_max) = 0, bound_tab_lev(nlev_max) = 0
   real(8) :: bound_tab_sp(nlev_max) = 0d0, bound_tab_cov(nlev_max) = 0d0

   !----------------------------------------------------------------------
   ! 正式选档规则（按初始候选粒子数 N、展开后最大轴向跨度 D 及当前模拟参数选择方法、网格和 spacing）
   !   * 10 <= N <= 1000   -> 连续软化直接求和（gid=0, 任务归 48 组）
   !   * 1000 < N <= 5000  -> 48/0.125, 64/0.125, 64/0.25, 128(0.125,0.25,0.5), 256 跨度表
   !   * 5000 < N <= 20000 -> 128/0.0625, 128/0.125；D>7.625 直接进入 256 序列
   !   * N > 20000         -> 256 跨度表（0.0625..4.0 + 严格 15 Mpc 终止档）
   !----------------------------------------------------------------------
   integer,parameter :: BOUND_SEL_OK = 0, BOUND_SEL_TOO_FEW = 1, BOUND_SEL_NO_COMBO = 2
   integer,parameter :: BOUND_METHOD_DIRECT = 0
   integer,parameter :: bound_fgl_nband = 3, bound_fgl_maxcand = 12, bound_fgl_maxlev = 32
   real(8),parameter :: bound_fgl_term_sentinel = -1d0
   integer :: bound_fgl_ncand(bound_fgl_nband) = 0
   integer :: bound_fgl_cgid(bound_fgl_nband, bound_fgl_maxcand) = 0
   real(8) :: bound_fgl_csp(bound_fgl_nband, bound_fgl_maxcand) = 0d0
   integer :: bound_fgl_nlev(4) = 0
   real(8) :: bound_fgl_sp(4, bound_fgl_maxlev) = 0d0
   real(8) :: bound_fgl_cov(4, bound_fgl_maxlev) = 0d0
   logical :: bound_fgl_ready = .false.

   ! 接入时由调用方提供粒子相空间、候选CRS与稳定ID。
   integer(8),allocatable :: pid(:)
   real(4),allocatable :: xv(:,:)
   integer,parameter :: bound_max_passes = 20
   integer,parameter :: bound_c2_min_np=20000, bound_c2_max_threads=16
   integer :: bound_c2_threads=1

   ! Continuous input is owned by the caller and never reordered by the solver.
   integer(4),allocatable :: bound_member_table(:),bound_member_offset(:),bound_member_count(:)

   ! Fixed, validated compact two-level configuration, not experimental defaults.
   integer,parameter :: bound_tl_nregions=512,bound_tl_n12=3926,bound_tl_n24=14688
   integer,parameter :: bound_tl_n_r2=128,bound_tl_n_r4=216,bound_green_threads=16
   integer,parameter :: bound_tl_max_slots=64
   real(8),parameter :: bound_tl_split=3.5d0
   integer :: bound_tl_slot_count=16
   type bound_inner_kernel_type
      integer :: n=0
      real(8) :: apm_soft=0d0,apm_range=0d0
      logical :: compact=.false.
      complex(4),allocatable :: green_phi(:,:,:)
      real(4) :: self_phi(-2:2,-2:2,-2:2)=0.0
   endtype
   type(bound_inner_kernel_type) :: bound_inner_gk_r2,bound_inner_gk_r4
   type bound_inner_slot_type
      integer(8) :: plan_r2=0,iplan_r2=0,plan_r4=0,iplan_r4=0
      real(4),pointer,contiguous :: rho_r2(:,:,:)=>null(),rho_r4(:,:,:)=>null()
      complex(4),pointer,contiguous :: ft_r2(:,:,:)=>null(),ft_r4(:,:,:)=>null()
   endtype
   type(bound_inner_slot_type) :: bound_inner_slots(bound_tl_max_slots)
   type bound_two_level_slot_ctx
      integer :: n_active=0
      integer :: active_regs(bound_tl_nregions)=0,region_r(bound_tl_nregions)=1
      integer,allocatable :: particle_r(:),target_reg(:),target_start(:),target_particles(:)
      integer,allocatable :: source_start(:),source_particles(:)
      real(4),allocatable :: self_in(:),wt_tgt(:,:,:),wt_src(:,:,:)
      integer,allocatable :: first_tsc_tgt(:,:),first_tsc_src(:,:)
   endtype
   type(bound_two_level_slot_ctx) :: bound_tl_ctx(fofcore)
   type bound_two_level_pp_type
      integer :: ncell(3)=0,ncell_total=0,max_search_k=1,ntasks=0
      real(8) :: u_min(3)=0d0,cw=0d0,max_cutoff=0d0
      integer,allocatable :: cell_id(:),cell_start(:),cell_particles(:)
      real(8),allocatable :: cell_max_a(:),task_max_a(:)
      integer,allocatable :: task_c1(:),task_p1(:),task_p2(:)
      real(4),allocatable :: phi_team(:,:)
   endtype
   type(bound_two_level_pp_type) :: bound_two_level_pps(fofcore)

contains

   subroutine bound_two_level_prepare_schedule(tasks)
      type(bound_schedule_type),intent(in) :: tasks
      integer :: h
      logical :: need_r2,need_r4
      need_r2=.false.; need_r4=.false.
      do h=1,size(tasks%method_gid)
         if (tasks%method_gid(h)/=4) cycle
         need_r2=need_r2.or.(bound_tl_split*tasks%spacing(h)/2d0>=dble(app))
         need_r4=need_r4.or.(bound_tl_split*tasks%spacing(h)/4d0>=dble(app))
      enddo
      if (need_r2.or.need_r4) then
         call bound_init_inner_gk(need_r2,need_r4)
         call bound_init_inner_slots(need_r2,need_r4)
      endif
   endsubroutine bound_two_level_prepare_schedule

   subroutine bound_init_inner_slots(need_r2,need_r4)
      logical,intent(in) :: need_r2,need_r4
      integer :: k,ios
      if (omp_in_parallel()) error stop 'bound_init_inner_slots: planner must be serial'
      if (bound_tl_slot_count<1 .or. bound_tl_slot_count>bound_tl_max_slots) &
         error stop 'bound_init_inner_slots: invalid slot count'
      do k=1,bound_tl_slot_count
         if (need_r2 .and. .not.associated(bound_inner_slots(k)%rho_r2)) then
            allocate(bound_inner_slots(k)%rho_r2(bound_tl_n_r2+2,bound_tl_n_r2,bound_tl_n_r2),stat=ios)
            if (ios/=0) error stop 'bound_init_inner_slots: r2 allocation failed'
            call c_f_pointer(c_loc(bound_inner_slots(k)%rho_r2),bound_inner_slots(k)%ft_r2, &
               [bound_tl_n_r2/2+1,bound_tl_n_r2,bound_tl_n_r2])
            call sfftw_plan_dft_r2c_3d(bound_inner_slots(k)%plan_r2,bound_tl_n_r2,bound_tl_n_r2,bound_tl_n_r2, &
               bound_inner_slots(k)%rho_r2,bound_inner_slots(k)%ft_r2,64)
            call sfftw_plan_dft_c2r_3d(bound_inner_slots(k)%iplan_r2,bound_tl_n_r2,bound_tl_n_r2,bound_tl_n_r2, &
               bound_inner_slots(k)%ft_r2,bound_inner_slots(k)%rho_r2,64)
            if (bound_inner_slots(k)%plan_r2==0 .or. bound_inner_slots(k)%iplan_r2==0) &
               error stop 'bound_init_inner_slots: r2 FFT plan failed'
         endif
         if (need_r4 .and. .not.associated(bound_inner_slots(k)%rho_r4)) then
            allocate(bound_inner_slots(k)%rho_r4(bound_tl_n_r4+2,bound_tl_n_r4,bound_tl_n_r4),stat=ios)
            if (ios/=0) error stop 'bound_init_inner_slots: r4 allocation failed'
            call c_f_pointer(c_loc(bound_inner_slots(k)%rho_r4),bound_inner_slots(k)%ft_r4, &
               [bound_tl_n_r4/2+1,bound_tl_n_r4,bound_tl_n_r4])
            call sfftw_plan_dft_r2c_3d(bound_inner_slots(k)%plan_r4,bound_tl_n_r4,bound_tl_n_r4,bound_tl_n_r4, &
               bound_inner_slots(k)%rho_r4,bound_inner_slots(k)%ft_r4,64)
            call sfftw_plan_dft_c2r_3d(bound_inner_slots(k)%iplan_r4,bound_tl_n_r4,bound_tl_n_r4,bound_tl_n_r4, &
               bound_inner_slots(k)%ft_r4,bound_inner_slots(k)%rho_r4,64)
            if (bound_inner_slots(k)%plan_r4==0 .or. bound_inner_slots(k)%iplan_r4==0) &
               error stop 'bound_init_inner_slots: r4 FFT plan failed'
         endif
      enddo
   endsubroutine bound_init_inner_slots

   subroutine bound_release_inner_slots
      integer :: k
      if (omp_in_parallel()) error stop 'bound_release_inner_slots: planner must be serial'
      do k=1,bound_tl_max_slots
         if (bound_inner_slots(k)%plan_r2/=0) call sfftw_destroy_plan(bound_inner_slots(k)%plan_r2)
         if (bound_inner_slots(k)%iplan_r2/=0) call sfftw_destroy_plan(bound_inner_slots(k)%iplan_r2)
         if (bound_inner_slots(k)%plan_r4/=0) call sfftw_destroy_plan(bound_inner_slots(k)%plan_r4)
         if (bound_inner_slots(k)%iplan_r4/=0) call sfftw_destroy_plan(bound_inner_slots(k)%iplan_r4)
         bound_inner_slots(k)%plan_r2=0; bound_inner_slots(k)%iplan_r2=0
         bound_inner_slots(k)%plan_r4=0; bound_inner_slots(k)%iplan_r4=0
         nullify(bound_inner_slots(k)%ft_r2,bound_inner_slots(k)%ft_r4)
         if (associated(bound_inner_slots(k)%rho_r2)) deallocate(bound_inner_slots(k)%rho_r2)
         if (associated(bound_inner_slots(k)%rho_r4)) deallocate(bound_inner_slots(k)%rho_r4)
      enddo
      bound_tl_slot_count=16
   endsubroutine bound_release_inner_slots

   subroutine bound_init_inner_gk(need_r2,need_r4)
      logical,intent(in) :: need_r2,need_r4
      integer :: ios
      if (allocated(bound_inner_gk_r2%green_phi)) then
         if (bound_inner_gk_r2%n/=bound_tl_n_r2 .or. bound_inner_gk_r2%apm_soft/=bound_tl_split .or. &
             bound_inner_gk_r2%apm_range/=7d0 .or. .not.bound_inner_gk_r2%compact) &
            deallocate(bound_inner_gk_r2%green_phi)
      endif
      if (allocated(bound_inner_gk_r4%green_phi)) then
         if (bound_inner_gk_r4%n/=bound_tl_n_r4 .or. bound_inner_gk_r4%apm_soft/=bound_tl_split .or. &
             bound_inner_gk_r4%apm_range/=14d0 .or. .not.bound_inner_gk_r4%compact) &
            deallocate(bound_inner_gk_r4%green_phi)
      endif
      if (need_r2 .and. .not.allocated(bound_inner_gk_r2%green_phi)) then
         call bound_init_inner_slots(.true.,.false.)
         allocate(bound_inner_gk_r2%green_phi(bound_tl_n_r2/2+1,bound_tl_n_r2,bound_tl_n_r2),stat=ios)
         if (ios/=0) error stop 'bound_init_inner_gk: r2 allocation failed'
         call bound_calc_green_and_self(bound_tl_n_r2,bound_tl_split,7d0,bound_inner_slots(1)%iplan_r2, &
            bound_inner_slots(1)%rho_r2,bound_inner_slots(1)%ft_r2, &
            bound_inner_gk_r2%green_phi,bound_inner_gk_r2%self_phi,truncated=.false.)
         bound_inner_gk_r2%n=bound_tl_n_r2
         bound_inner_gk_r2%apm_soft=bound_tl_split; bound_inner_gk_r2%apm_range=7d0
         bound_inner_gk_r2%compact=.true.
      endif
      if (need_r4 .and. .not.allocated(bound_inner_gk_r4%green_phi)) then
         call bound_init_inner_slots(.false.,.true.)
         allocate(bound_inner_gk_r4%green_phi(bound_tl_n_r4/2+1,bound_tl_n_r4,bound_tl_n_r4),stat=ios)
         if (ios/=0) error stop 'bound_init_inner_gk: r4 allocation failed'
         call bound_calc_green_and_self(bound_tl_n_r4,bound_tl_split,14d0,bound_inner_slots(1)%iplan_r4, &
            bound_inner_slots(1)%rho_r4,bound_inner_slots(1)%ft_r4, &
            bound_inner_gk_r4%green_phi,bound_inner_gk_r4%self_phi,truncated=.false.)
         bound_inner_gk_r4%n=bound_tl_n_r4
         bound_inner_gk_r4%apm_soft=bound_tl_split; bound_inner_gk_r4%apm_range=14d0
         bound_inner_gk_r4%compact=.true.
      endif
   endsubroutine bound_init_inner_gk

   pure subroutine bound_region_bounds(b,lo,hi)
      integer,intent(in) :: b
      real(4),intent(out) :: lo(3),hi(3)
      integer :: bx,by,bz
      bx=mod(b-1,8); by=mod((b-1)/8,8); bz=(b-1)/64
      lo(1)=real(bx*16,4); hi(1)=min(lo(1)+16.0_4,127.0_4)
      lo(2)=real(by*16,4); hi(2)=min(lo(2)+16.0_4,127.0_4)
      lo(3)=real(bz*16,4); hi(3)=min(lo(3)+16.0_4,127.0_4)
   endsubroutine bound_region_bounds

   pure integer function bound_get_region_id(coord_grid) result(rid)
      real(4),intent(in) :: coord_grid(3)
      integer :: bx,by,bz,cell(3)
      cell=floor(dble(coord_grid))
      bx=min(7,max(0,cell(1)/16))
      by=min(7,max(0,cell(2)/16))
      bz=min(7,max(0,cell(3)/16))
      rid=1+bx+by*8+bz*64
   endfunction bound_get_region_id

   subroutine bound_two_level_setup_ratios(spacing,app_val,region_r,particle_r,target_reg,q,np,origin)
      integer,intent(in) :: np
      real(4),intent(in) :: spacing,app_val,origin(3),q(3,np)
      integer,intent(out) :: region_r(bound_tl_nregions),particle_r(np),target_reg(np)
      integer :: b,i,r_req,r_act,b_cnt(bound_tl_nregions)
      real(4) :: t(3)
      real(8) :: h_b
      b_cnt=0; region_r=1
      do i=1,np
         t=real((dble(q(:,i))-dble(origin))/dble(spacing),4)
         b=bound_get_region_id(t)
         if (b<1 .or. b>bound_tl_nregions) error stop 'two-level: target region outside configuration'
         target_reg(i)=b
         b_cnt(b)=b_cnt(b)+1
      enddo
      do b=1,bound_tl_nregions
         if (b_cnt(b)<bound_tl_n12) then
            r_req=1
         else if (b_cnt(b)<bound_tl_n24) then
            r_req=2
         else
            r_req=4
         endif
         r_act=r_req
         if (r_act==4) then
            h_b=dble(spacing)/4d0
            if (bound_tl_split*h_b<dble(app_val)) r_act=2
         endif
         if (r_act==2) then
            h_b=dble(spacing)/2d0
            if (bound_tl_split*h_b<dble(app_val)) r_act=1
         endif
         region_r(b)=r_act
      enddo
      do i=1,np
         particle_r(i)=region_r(target_reg(i))
      enddo
   endsubroutine bound_two_level_setup_ratios

   subroutine bound_two_level_release_halo(s)
      integer,intent(in) :: s
      if (s<1 .or. s>fofcore) error stop 'bound_two_level_release_halo: invalid owner'
      if (allocated(bound_tl_ctx(s)%particle_r)) deallocate(bound_tl_ctx(s)%particle_r)
      if (allocated(bound_tl_ctx(s)%target_reg)) deallocate(bound_tl_ctx(s)%target_reg)
      if (allocated(bound_tl_ctx(s)%target_start)) deallocate(bound_tl_ctx(s)%target_start)
      if (allocated(bound_tl_ctx(s)%target_particles)) deallocate(bound_tl_ctx(s)%target_particles)
      if (allocated(bound_tl_ctx(s)%source_start)) deallocate(bound_tl_ctx(s)%source_start)
      if (allocated(bound_tl_ctx(s)%source_particles)) deallocate(bound_tl_ctx(s)%source_particles)
      if (allocated(bound_tl_ctx(s)%self_in)) deallocate(bound_tl_ctx(s)%self_in)
      if (allocated(bound_tl_ctx(s)%wt_tgt)) deallocate(bound_tl_ctx(s)%wt_tgt)
      if (allocated(bound_tl_ctx(s)%first_tsc_tgt)) deallocate(bound_tl_ctx(s)%first_tsc_tgt)
      if (allocated(bound_tl_ctx(s)%wt_src)) deallocate(bound_tl_ctx(s)%wt_src)
      if (allocated(bound_tl_ctx(s)%first_tsc_src)) deallocate(bound_tl_ctx(s)%first_tsc_src)
      bound_tl_ctx(s)%n_active=0
      bound_tl_ctx(s)%region_r=1
   endsubroutine bound_two_level_release_halo

   subroutine bound_two_level_build_halo(s,np,q,spacing,origin,app_val)
      integer,intent(in) :: s,np
      real(4),intent(in) :: q(3,np),spacing,origin(3),app_val
      integer :: b,r_b,n_in,i,j,axis,ios,n_act,act_reg,p,total_src
      integer :: target_cnt(bound_tl_nregions),source_cnt(bound_tl_nregions)
      real(4) :: u_min(3),u_max(3),b_outer,u0_reg(3),t_fine(3),delta,coord,sw,tw
      real(4),allocatable :: t(:,:)
      integer :: nearest,sz,sy,sx,tz,ty,tx,src_lo(3),src_hi(3),tgt_lo(3),tgt_hi(3),axis_gap(3)
      real(8) :: acc_self
      logical :: need_r2,need_r4
      if (s<1 .or. s>fofcore) error stop 'bound_two_level_build_halo: invalid owner'
      call bound_two_level_release_halo(s)
      allocate(bound_tl_ctx(s)%particle_r(np),bound_tl_ctx(s)%target_reg(np),stat=ios)
      if (ios/=0) error stop 'bound_two_level_build_halo: ratio allocation failed'
      call bound_two_level_setup_ratios(spacing,app_val,bound_tl_ctx(s)%region_r, &
         bound_tl_ctx(s)%particle_r,bound_tl_ctx(s)%target_reg,q,np,origin)
      target_cnt=0
      do i=1,np
         b=bound_tl_ctx(s)%target_reg(i)
         target_cnt(b)=target_cnt(b)+1
      enddo
      if (all(bound_tl_ctx(s)%particle_r==1)) return
      n_act=0
      do b=1,bound_tl_nregions
         if (target_cnt(b)>0 .and. bound_tl_ctx(s)%region_r(b)>1) then
            n_act=n_act+1
            bound_tl_ctx(s)%active_regs(n_act)=b
         endif
      enddo
      bound_tl_ctx(s)%n_active=n_act
      if (n_act==0) return
      need_r2=any(bound_tl_ctx(s)%particle_r==2)
      need_r4=any(bound_tl_ctx(s)%particle_r==4)
      allocate(t(3,np),stat=ios)
      if (ios/=0) error stop 'bound_two_level_build_halo: coordinate allocation failed'
      if (need_r2 .and. .not.allocated(bound_inner_gk_r2%green_phi)) error stop 'r2 kernel not prepared'
      if (need_r4 .and. .not.allocated(bound_inner_gk_r4%green_phi)) error stop 'r4 kernel not prepared'
      allocate(bound_tl_ctx(s)%target_start(bound_tl_nregions+1),bound_tl_ctx(s)%target_particles(np),stat=ios)
      if (ios/=0) error stop 'bound_two_level_build_halo: target CSR allocation failed'
      bound_tl_ctx(s)%target_start(1)=1
      do b=1,bound_tl_nregions
         bound_tl_ctx(s)%target_start(b+1)=bound_tl_ctx(s)%target_start(b)+target_cnt(b)
      enddo
      target_cnt=0
      do i=1,np
         b=bound_tl_ctx(s)%target_reg(i)
         p=bound_tl_ctx(s)%target_start(b)+target_cnt(b)
         bound_tl_ctx(s)%target_particles(p)=i
         target_cnt(b)=target_cnt(b)+1
      enddo
      do i=1,np
         t(:,i)=real((dble(q(:,i))-dble(origin))/dble(spacing),4)
      enddo
      source_cnt=0
      do act_reg=1,n_act
         b=bound_tl_ctx(s)%active_regs(act_reg)
         r_b=bound_tl_ctx(s)%region_r(b)
         call bound_region_bounds(b,u_min,u_max)
         b_outer=3.5+3.0/real(r_b,4)
         do j=1,np
            if (t(1,j)>=u_min(1)-b_outer .and. t(1,j)<=u_max(1)+b_outer .and. &
                t(2,j)>=u_min(2)-b_outer .and. t(2,j)<=u_max(2)+b_outer .and. &
                t(3,j)>=u_min(3)-b_outer .and. t(3,j)<=u_max(3)+b_outer) source_cnt(b)=source_cnt(b)+1
         enddo
      enddo
      allocate(bound_tl_ctx(s)%source_start(bound_tl_nregions+1),stat=ios)
      if (ios/=0) error stop 'bound_two_level_build_halo: source offsets allocation failed'
      bound_tl_ctx(s)%source_start(1)=1
      do b=1,bound_tl_nregions
         bound_tl_ctx(s)%source_start(b+1)=bound_tl_ctx(s)%source_start(b)+source_cnt(b)
      enddo
      total_src=bound_tl_ctx(s)%source_start(bound_tl_nregions+1)-1
      allocate(bound_tl_ctx(s)%source_particles(max(1,total_src)), &
         bound_tl_ctx(s)%wt_src(0:2,3,max(1,total_src)), &
         bound_tl_ctx(s)%first_tsc_src(3,max(1,total_src)),stat=ios)
      if (ios/=0) error stop 'bound_two_level_build_halo: source CSR allocation failed'
      source_cnt=0
      do act_reg=1,n_act
         b=bound_tl_ctx(s)%active_regs(act_reg)
         r_b=bound_tl_ctx(s)%region_r(b)
         n_in=merge(bound_tl_n_r2,bound_tl_n_r4,r_b==2)
         call bound_region_bounds(b,u_min,u_max)
         b_outer=3.5+3.0/real(r_b,4)
         u0_reg=u_min-b_outer-(1.0/real(r_b,4))+(0.0/real(r_b,4))
         do j=1,np
            if (t(1,j)>=u_min(1)-b_outer .and. t(1,j)<=u_max(1)+b_outer .and. &
                t(2,j)>=u_min(2)-b_outer .and. t(2,j)<=u_max(2)+b_outer .and. &
                t(3,j)>=u_min(3)-b_outer .and. t(3,j)<=u_max(3)+b_outer) then
               p=bound_tl_ctx(s)%source_start(b)+source_cnt(b)
               bound_tl_ctx(s)%source_particles(p)=j
               t_fine=(t(:,j)-u0_reg)*real(r_b,4)
               do axis=1,3
                  coord=t_fine(axis)
                  nearest=nint(coord)
                  delta=coord-real(nearest,4)
                  bound_tl_ctx(s)%first_tsc_src(axis,p)=nearest
                  bound_tl_ctx(s)%wt_src(0,axis,p)=0.5*(0.5-delta)**2
                  bound_tl_ctx(s)%wt_src(1,axis,p)=0.75-delta**2
                  bound_tl_ctx(s)%wt_src(2,axis,p)=0.5*(0.5+delta)**2
                  if (nearest<1 .or. nearest+2>n_in/2) error stop 'source node out of half domain'
               enddo
               source_cnt(b)=source_cnt(b)+1
            endif
         enddo
      enddo
      allocate(bound_tl_ctx(s)%self_in(np),bound_tl_ctx(s)%wt_tgt(0:2,3,np), &
         bound_tl_ctx(s)%first_tsc_tgt(3,np),stat=ios)
      if (ios/=0) error stop 'bound_two_level_build_halo: target arrays allocation failed'
      bound_tl_ctx(s)%self_in=0.0
      do act_reg=1,n_act
         b=bound_tl_ctx(s)%active_regs(act_reg)
         r_b=bound_tl_ctx(s)%region_r(b)
         n_in=merge(bound_tl_n_r2,bound_tl_n_r4,r_b==2)
         call bound_region_bounds(b,u_min,u_max)
         b_outer=3.5+3.0/real(r_b,4)
         u0_reg=u_min-b_outer-(1.0/real(r_b,4))+(0.0/real(r_b,4))
         do p=bound_tl_ctx(s)%target_start(b),bound_tl_ctx(s)%target_start(b+1)-1
            i=bound_tl_ctx(s)%target_particles(p)
            t_fine=(t(:,i)-u0_reg)*real(r_b,4)
            do axis=1,3
               coord=t_fine(axis)
               nearest=nint(coord)
               delta=coord-real(nearest,4)
               bound_tl_ctx(s)%first_tsc_tgt(axis,i)=nearest
               bound_tl_ctx(s)%wt_tgt(0,axis,i)=0.5*(0.5-delta)**2
               bound_tl_ctx(s)%wt_tgt(1,axis,i)=0.75-delta**2
               bound_tl_ctx(s)%wt_tgt(2,axis,i)=0.5*(0.5+delta)**2
               if (nearest<1 .or. nearest+2>n_in/2) error stop 'target node out of half domain'
            enddo
            acc_self=0d0
            do sz=0,2
               do sy=0,2
                  do sx=0,2
                     sw=bound_tl_ctx(s)%wt_tgt(sx,1,i)*bound_tl_ctx(s)%wt_tgt(sy,2,i)* &
                        bound_tl_ctx(s)%wt_tgt(sz,3,i)
                     do tz=0,2
                        do ty=0,2
                           do tx=0,2
                              tw=bound_tl_ctx(s)%wt_tgt(tx,1,i)*bound_tl_ctx(s)%wt_tgt(ty,2,i)* &
                                 bound_tl_ctx(s)%wt_tgt(tz,3,i)
                              if (r_b==2) then
                                 acc_self=acc_self+dble(sw*tw)*dble(bound_inner_gk_r2%self_phi(tx-sx,ty-sy,tz-sz))
                              else
                                 acc_self=acc_self+dble(sw*tw)*dble(bound_inner_gk_r4%self_phi(tx-sx,ty-sy,tz-sz))
                              endif
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            bound_tl_ctx(s)%self_in(i)=real(acc_self,4)
         enddo
      enddo
      do act_reg=1,n_act
         b=bound_tl_ctx(s)%active_regs(act_reg)
         r_b=bound_tl_ctx(s)%region_r(b)
         n_in=merge(bound_tl_n_r2,bound_tl_n_r4,r_b==2)
         src_lo=huge(1); src_hi=-huge(1)
         tgt_lo=huge(1); tgt_hi=-huge(1)
         do p=bound_tl_ctx(s)%source_start(b),bound_tl_ctx(s)%source_start(b+1)-1
            src_lo=min(src_lo,bound_tl_ctx(s)%first_tsc_src(:,p))
            src_hi=max(src_hi,bound_tl_ctx(s)%first_tsc_src(:,p)+2)
         enddo
         do p=bound_tl_ctx(s)%target_start(b),bound_tl_ctx(s)%target_start(b+1)-1
            i=bound_tl_ctx(s)%target_particles(p)
            tgt_lo=min(tgt_lo,bound_tl_ctx(s)%first_tsc_tgt(:,i))
            tgt_hi=max(tgt_hi,bound_tl_ctx(s)%first_tsc_tgt(:,i)+2)
         enddo
         axis_gap=max(abs(tgt_hi-src_lo),abs(src_hi-tgt_lo))
         if (any(src_lo<1) .or. any(src_hi>n_in/2) .or. any(tgt_lo<1) .or. any(tgt_hi>n_in/2) .or. &
             any(axis_gap>=n_in/2)) error stop 'two-level compact source/target nodes exceed first-half support'
      enddo
   endsubroutine bound_two_level_build_halo

   subroutine bound_two_level_run_pm(s,np,keep,spacing,a_val,phi_in,nth)
      integer,intent(in) :: s,np,nth
      logical,intent(in) :: keep(np)
      real(4),intent(in) :: spacing
      real(8),intent(in) :: a_val
      real(4),intent(out) :: phi_in(np)
      integer :: n_act,reg_idx,b,r_b,n_in,islot,tid,p,j,i,team_stride,sz,sy,sx,iz,iy,ix
      real(4) :: raw_in
      phi_in=0.0
      n_act=bound_tl_ctx(s)%n_active
      if (n_act==0) return
      team_stride=max(1,bound_c2_threads,nth)
      if (nth<1 .or. nth>4 .or. s<1 .or. s>fofcore .or. &
          (s-1)*team_stride+nth>bound_tl_slot_count) error stop 'two_level_run_pm: invalid owner/team'
      !$omp parallel num_threads(nth) default(shared) &
      !$omp private(tid,islot,reg_idx,b,r_b,n_in,p,j,i,sz,sy,sx,iz,iy,ix,raw_in)
      !$omp single
      if (omp_get_num_threads()/=nth) error stop 'two_level_run_pm: actual team size mismatch'
      !$omp end single
      tid=omp_get_thread_num()+1
      islot=(s-1)*team_stride+tid
      !$omp do schedule(dynamic)
      do reg_idx=1,n_act
         b=bound_tl_ctx(s)%active_regs(reg_idx)
         r_b=bound_tl_ctx(s)%region_r(b)
         n_in=merge(bound_tl_n_r2,bound_tl_n_r4,r_b==2)
         if (r_b==2) then
            bound_inner_slots(islot)%rho_r2=0.0
            do p=bound_tl_ctx(s)%source_start(b),bound_tl_ctx(s)%source_start(b+1)-1
               j=bound_tl_ctx(s)%source_particles(p)
               if (.not.keep(j)) cycle
               do sz=0,2
                  iz=bound_tl_ctx(s)%first_tsc_src(3,p)+sz
                  do sy=0,2
                     iy=bound_tl_ctx(s)%first_tsc_src(2,p)+sy
                     do sx=0,2
                        ix=bound_tl_ctx(s)%first_tsc_src(1,p)+sx
                        bound_inner_slots(islot)%rho_r2(ix,iy,iz)=bound_inner_slots(islot)%rho_r2(ix,iy,iz)+ &
                           bound_tl_ctx(s)%wt_src(sx,1,p)*bound_tl_ctx(s)%wt_src(sy,2,p)*bound_tl_ctx(s)%wt_src(sz,3,p)
                     enddo
                  enddo
               enddo
            enddo
         else
            bound_inner_slots(islot)%rho_r4=0.0
            do p=bound_tl_ctx(s)%source_start(b),bound_tl_ctx(s)%source_start(b+1)-1
               j=bound_tl_ctx(s)%source_particles(p)
               if (.not.keep(j)) cycle
               do sz=0,2
                  iz=bound_tl_ctx(s)%first_tsc_src(3,p)+sz
                  do sy=0,2
                     iy=bound_tl_ctx(s)%first_tsc_src(2,p)+sy
                     do sx=0,2
                        ix=bound_tl_ctx(s)%first_tsc_src(1,p)+sx
                        bound_inner_slots(islot)%rho_r4(ix,iy,iz)=bound_inner_slots(islot)%rho_r4(ix,iy,iz)+ &
                           bound_tl_ctx(s)%wt_src(sx,1,p)*bound_tl_ctx(s)%wt_src(sy,2,p)*bound_tl_ctx(s)%wt_src(sz,3,p)
                     enddo
                  enddo
               enddo
            enddo
         endif
         if (r_b==2) then
            call sfftw_execute_dft_r2c(bound_inner_slots(islot)%plan_r2, &
               bound_inner_slots(islot)%rho_r2,bound_inner_slots(islot)%ft_r2)
            bound_inner_slots(islot)%ft_r2=bound_inner_slots(islot)%ft_r2*bound_inner_gk_r2%green_phi
            call sfftw_execute_dft_c2r(bound_inner_slots(islot)%iplan_r2, &
               bound_inner_slots(islot)%ft_r2,bound_inner_slots(islot)%rho_r2)
            bound_inner_slots(islot)%rho_r2=bound_inner_slots(islot)%rho_r2/real(n_in**3,4)
         else
            call sfftw_execute_dft_r2c(bound_inner_slots(islot)%plan_r4, &
               bound_inner_slots(islot)%rho_r4,bound_inner_slots(islot)%ft_r4)
            bound_inner_slots(islot)%ft_r4=bound_inner_slots(islot)%ft_r4*bound_inner_gk_r4%green_phi
            call sfftw_execute_dft_c2r(bound_inner_slots(islot)%iplan_r4, &
               bound_inner_slots(islot)%ft_r4,bound_inner_slots(islot)%rho_r4)
            bound_inner_slots(islot)%rho_r4=bound_inner_slots(islot)%rho_r4/real(n_in**3,4)
         endif
         do p=bound_tl_ctx(s)%target_start(b),bound_tl_ctx(s)%target_start(b+1)-1
            i=bound_tl_ctx(s)%target_particles(p)
            if (.not.keep(i)) cycle
            raw_in=0.0
            if (r_b==2) then
               do sz=0,2
                  iz=bound_tl_ctx(s)%first_tsc_tgt(3,i)+sz
                  do sy=0,2
                     iy=bound_tl_ctx(s)%first_tsc_tgt(2,i)+sy
                     do sx=0,2
                        ix=bound_tl_ctx(s)%first_tsc_tgt(1,i)+sx
                        raw_in=raw_in+bound_inner_slots(islot)%rho_r2(ix,iy,iz)*bound_tl_ctx(s)%wt_tgt(sx,1,i)* &
                           bound_tl_ctx(s)%wt_tgt(sy,2,i)*bound_tl_ctx(s)%wt_tgt(sz,3,i)
                     enddo
                  enddo
               enddo
            else
               do sz=0,2
                  iz=bound_tl_ctx(s)%first_tsc_tgt(3,i)+sz
                  do sy=0,2
                     iy=bound_tl_ctx(s)%first_tsc_tgt(2,i)+sy
                     do sx=0,2
                        ix=bound_tl_ctx(s)%first_tsc_tgt(1,i)+sx
                        raw_in=raw_in+bound_inner_slots(islot)%rho_r4(ix,iy,iz)*bound_tl_ctx(s)%wt_tgt(sx,1,i)* &
                           bound_tl_ctx(s)%wt_tgt(sy,2,i)*bound_tl_ctx(s)%wt_tgt(sz,3,i)
                     enddo
                  enddo
               enddo
            endif
            phi_in(i)=real((a_val*dble(r_b)/(6d0*pi*dble(spacing))),4)*(raw_in-bound_tl_ctx(s)%self_in(i))
         enddo
      enddo
      !$omp end do
      !$omp end parallel
   endsubroutine bound_two_level_run_pm

   subroutine bound_two_level_pp_build(s,np,q,spacing,ratios)
      integer,intent(in) :: s,np,ratios(np)
      real(4),intent(in) :: q(3,np),spacing
      integer :: i,c,nt,ic(3),ncell(3),need,fill_ptr,ntask,p_cnt,p_cur,p_next
      integer :: r_max,r_min,r_val,k_search
      real(8) :: u_min(3),u_max(3),cw,a_val_i,a_max_halo,a_max_task,pos(3),a_table(4)
      if (s<1 .or. s>fofcore .or. np<=0) error stop 'two_level_pp_build: invalid slot or count'
      if (.not.ieee_is_finite(spacing) .or. spacing<=0.0) error stop 'two_level_pp_build: invalid spacing'
      if (.not.all(ieee_is_finite(q))) error stop 'two_level_pp_build: nonfinite position'
      if (any(ratios/=1 .and. ratios/=2 .and. ratios/=4)) error stop 'two_level_pp_build: invalid ratios'
      r_max=1; r_min=4
      do i=1,np
         r_val=ratios(i)
         if (r_val<1) r_val=1
         if (r_val>4) r_val=4
         if (r_val>r_max) r_max=r_val
         if (r_val<r_min) r_min=r_val
      enddo
      a_table(1)=3.5d0*dble(spacing)
      a_table(2)=bound_tl_split*dble(spacing)/2d0
      a_table(3)=0d0
      a_table(4)=bound_tl_split*dble(spacing)/4d0
      if (r_max==1) then
         cw=3.5d0*dble(spacing)
      else
         cw=bound_tl_split*dble(spacing)/dble(r_max)
      endif
      if (cw<=0d0) cw=3.5d0*dble(spacing)
      if (r_min==1) then
         a_max_halo=3.5d0*dble(spacing)
      else
         a_max_halo=bound_tl_split*dble(spacing)/dble(r_min)
      endif
      if (any(ratios==1)) a_max_halo=max(a_max_halo,3.5d0*dble(spacing))
      k_search=max(1,min(4,int(ceiling(a_max_halo/cw))))
      u_min=1d30; u_max=-1d30
      do i=1,np
         pos=dble(q(:,i))
         u_min=min(u_min,pos); u_max=max(u_max,pos)
      enddo
      do i=1,3
         ncell(i)=max(1,int(floor((u_max(i)-u_min(i))/cw))+1)
      enddo
      nt=ncell(1)*ncell(2)*ncell(3)
      need=max(np,1)
      if (.not.allocated(bound_two_level_pps(s)%cell_id)) then
         allocate(bound_two_level_pps(s)%cell_id(need),bound_two_level_pps(s)%cell_particles(need))
      else if (size(bound_two_level_pps(s)%cell_id)<need) then
         deallocate(bound_two_level_pps(s)%cell_id,bound_two_level_pps(s)%cell_particles)
         allocate(bound_two_level_pps(s)%cell_id(need),bound_two_level_pps(s)%cell_particles(need))
      endif
      if (.not.allocated(bound_two_level_pps(s)%cell_start)) then
         allocate(bound_two_level_pps(s)%cell_start(nt+1),bound_two_level_pps(s)%cell_max_a(nt))
      else if (size(bound_two_level_pps(s)%cell_start)<nt+1) then
         deallocate(bound_two_level_pps(s)%cell_start,bound_two_level_pps(s)%cell_max_a)
         allocate(bound_two_level_pps(s)%cell_start(nt+1),bound_two_level_pps(s)%cell_max_a(nt))
      endif
      bound_two_level_pps(s)%ncell=ncell
      bound_two_level_pps(s)%ncell_total=nt
      bound_two_level_pps(s)%u_min=u_min
      bound_two_level_pps(s)%cw=cw
      bound_two_level_pps(s)%max_search_k=k_search
      bound_two_level_pps(s)%max_cutoff=a_max_halo
      bound_two_level_pps(s)%cell_start(1:nt+1)=0
      bound_two_level_pps(s)%cell_max_a(1:nt)=0d0
      do i=1,np
         pos=dble(q(:,i))
         ic=min(max(int(floor((pos-u_min)/cw))+1,1),ncell)
         c=ic(1)+(ic(2)-1)*ncell(1)+(ic(3)-1)*ncell(1)*ncell(2)
         bound_two_level_pps(s)%cell_id(i)=c
         bound_two_level_pps(s)%cell_start(c+1)=bound_two_level_pps(s)%cell_start(c+1)+1
         r_val=max(1,min(4,ratios(i)))
         a_val_i=a_table(r_val)
         if (a_val_i>bound_two_level_pps(s)%cell_max_a(c)) bound_two_level_pps(s)%cell_max_a(c)=a_val_i
      enddo
      do c=1,nt
         bound_two_level_pps(s)%cell_start(c+1)= &
            bound_two_level_pps(s)%cell_start(c+1)+bound_two_level_pps(s)%cell_start(c)
      enddo
      bound_two_level_pps(s)%cell_start(1)=1
      do c=2,nt+1
         bound_two_level_pps(s)%cell_start(c)=bound_two_level_pps(s)%cell_start(c)+1
      enddo
      do i=1,np
         c=bound_two_level_pps(s)%cell_id(i)
         fill_ptr=bound_two_level_pps(s)%cell_start(c)
         bound_two_level_pps(s)%cell_particles(fill_ptr)=i
         bound_two_level_pps(s)%cell_start(c)=fill_ptr+1
      enddo
      do c=nt,2,-1
         bound_two_level_pps(s)%cell_start(c)=bound_two_level_pps(s)%cell_start(c-1)
      enddo
      bound_two_level_pps(s)%cell_start(1)=1
      ntask=0
      do c=1,nt
         p_cnt=bound_two_level_pps(s)%cell_start(c+1)-bound_two_level_pps(s)%cell_start(c)
         if (p_cnt<=0) cycle
         ntask=ntask+(p_cnt+1023)/1024
      enddo
      if (.not.allocated(bound_two_level_pps(s)%task_c1)) then
         allocate(bound_two_level_pps(s)%task_c1(max(ntask,1)),bound_two_level_pps(s)%task_p1(max(ntask,1)), &
            bound_two_level_pps(s)%task_p2(max(ntask,1)),bound_two_level_pps(s)%task_max_a(max(ntask,1)))
      else if (size(bound_two_level_pps(s)%task_c1)<ntask) then
         deallocate(bound_two_level_pps(s)%task_c1,bound_two_level_pps(s)%task_p1, &
            bound_two_level_pps(s)%task_p2,bound_two_level_pps(s)%task_max_a)
         allocate(bound_two_level_pps(s)%task_c1(ntask),bound_two_level_pps(s)%task_p1(ntask), &
            bound_two_level_pps(s)%task_p2(ntask),bound_two_level_pps(s)%task_max_a(ntask))
      endif
      bound_two_level_pps(s)%ntasks=ntask
      ntask=0
      do c=1,nt
         p_cnt=bound_two_level_pps(s)%cell_start(c+1)-bound_two_level_pps(s)%cell_start(c)
         if (p_cnt<=0) cycle
         p_cur=bound_two_level_pps(s)%cell_start(c)
         do while (p_cur<=bound_two_level_pps(s)%cell_start(c+1)-1)
            p_next=min(p_cur+1023,bound_two_level_pps(s)%cell_start(c+1)-1)
            ntask=ntask+1
            bound_two_level_pps(s)%task_c1(ntask)=c
            bound_two_level_pps(s)%task_p1(ntask)=p_cur
            bound_two_level_pps(s)%task_p2(ntask)=p_next
            a_max_task=0d0
            do i=p_cur,p_next
               r_val=max(1,min(4,ratios(bound_two_level_pps(s)%cell_particles(i))))
               a_val_i=a_table(r_val)
               if (a_val_i>a_max_task) a_max_task=a_val_i
            enddo
            bound_two_level_pps(s)%task_max_a(ntask)=a_max_task
            p_cur=p_next+1
         enddo
      enddo
   endsubroutine bound_two_level_pp_build

   subroutine bound_two_level_pp_full(s,np,q,keep,spacing,app_val,a_val,ratios,phi_c2,nth)
      integer,intent(in) :: s,np,nth,ratios(np)
      real(4),intent(in) :: q(3,np),spacing,app_val
      logical,intent(in) :: keep(np)
      real(8),intent(in) :: a_val
      real(4),intent(inout) :: phi_c2(np)
      integer :: n_k,i,p1,p2,i1,i2,c1,c2,ntasks,itask,p1_beg,p1_end
      integer :: cx1,cy1,cz1,dx,dy,dz,nx,ny,nz,gap_x,gap_y,gap_z,nc(3),k_search,tid,tid_k
      real(8) :: cw,cw2,dmin2,rcut,task_max_a,eps,factor,dx_r,dy_r,dz_r,r2,d,phi_eps
      real(8) :: a_table(4),a1,a2,r2_max_12,delta_depth1,delta_depth2
      real(4) :: dphi1,dphi2,phi_sum
      if (nth<1 .or. nth>4) error stop 'two_level_pp_full: invalid team size'
      if (s<1 .or. s>fofcore .or. np<=0) error stop 'two_level_pp_full: invalid slot or count'
      if (spacing<=0.0 .or. app_val<=0.0) error stop 'two_level_pp_full: invalid length'
      phi_c2=0.0
      n_k=count(keep(1:np))
      if (n_k<=1) return
      a_table(1)=3.5d0*dble(spacing)/1d0
      a_table(2)=bound_tl_split*dble(spacing)/2d0
      a_table(3)=0d0
      a_table(4)=bound_tl_split*dble(spacing)/4d0
      eps=dble(app_val)
      factor=a_val/(6d0*pi)
      nc=bound_two_level_pps(s)%ncell
      cw=bound_two_level_pps(s)%cw
      cw2=cw*cw
      k_search=bound_two_level_pps(s)%max_search_k
      ntasks=bound_two_level_pps(s)%ntasks
      if (nth>1) then
         if (.not.allocated(bound_two_level_pps(s)%phi_team)) then
            allocate(bound_two_level_pps(s)%phi_team(max(np,1),nth))
         else if (size(bound_two_level_pps(s)%phi_team,1)<np .or. &
                  size(bound_two_level_pps(s)%phi_team,2)<nth) then
            deallocate(bound_two_level_pps(s)%phi_team)
            allocate(bound_two_level_pps(s)%phi_team(max(np,1),nth))
         endif
         bound_two_level_pps(s)%phi_team(1:np,1:nth)=0.0
         !$omp parallel num_threads(nth) default(shared) &
         !$omp private(tid,tid_k,itask,c1,cx1,cy1,cz1,p1_beg,p1_end,p1,p2,i1,i2, &
         !$omp dx_r,dy_r,dz_r,r2,d,phi_eps,dmin2,rcut,task_max_a,gap_x,gap_y,gap_z, &
         !$omp dz,dy,dx,nx,ny,nz,c2,delta_depth1,delta_depth2,dphi1,dphi2,a1,a2,r2_max_12,i,phi_sum)
         !$omp single
         if (omp_get_num_threads()/=nth) error stop 'two_level_pp_full: actual team size mismatch'
         !$omp end single
         tid=omp_get_thread_num()+1
         !$omp do schedule(dynamic,1)
         do itask=1,ntasks
            c1=bound_two_level_pps(s)%task_c1(itask)
            p1_beg=bound_two_level_pps(s)%task_p1(itask)
            p1_end=bound_two_level_pps(s)%task_p2(itask)
            task_max_a=bound_two_level_pps(s)%task_max_a(itask)
            cz1=(c1-1)/(nc(1)*nc(2))+1
            cy1=mod(c1-1,nc(1)*nc(2))/nc(1)+1
            cx1=mod(mod(c1-1,nc(1)*nc(2)),nc(1))+1
            do p1=p1_beg,p1_end
               i1=bound_two_level_pps(s)%cell_particles(p1)
               if (.not.keep(i1)) cycle
               a1=a_table(max(1,min(4,ratios(i1))))
               do p2=p1+1,bound_two_level_pps(s)%cell_start(c1+1)-1
                  i2=bound_two_level_pps(s)%cell_particles(p2)
                  if (.not.keep(i2)) cycle
                  a2=a_table(max(1,min(4,ratios(i2))))
                  dx_r=dble(q(1,i2))-dble(q(1,i1))
                  dy_r=dble(q(2,i2))-dble(q(2,i1))
                  dz_r=dble(q(3,i2))-dble(q(3,i1))
                  r2=dx_r*dx_r+dy_r*dy_r+dz_r*dz_r
                  r2_max_12=max(a1,a2)*max(a1,a2)
                  if (r2<r2_max_12) then
                     d=sqrt(r2)
                     phi_eps=bound_depth_continuous(d,eps)
                     if (r2<a1*a1) then
                        delta_depth1=phi_eps-bound_depth_continuous(d,a1)
                        dphi1=-real(factor*delta_depth1,4)
                        bound_two_level_pps(s)%phi_team(i1,tid)=bound_two_level_pps(s)%phi_team(i1,tid)+dphi1
                     endif
                     if (r2<a2*a2) then
                        delta_depth2=phi_eps-bound_depth_continuous(d,a2)
                        dphi2=-real(factor*delta_depth2,4)
                        bound_two_level_pps(s)%phi_team(i2,tid)=bound_two_level_pps(s)%phi_team(i2,tid)+dphi2
                     endif
                  endif
               enddo
            enddo
            do dz=0,k_search
               nz=cz1+dz
               if (nz<1 .or. nz>nc(3)) cycle
               gap_z=max(0,abs(dz)-1)
               do dy=-k_search,k_search
                  ny=cy1+dy
                  if (ny<1 .or. ny>nc(2)) cycle
                  gap_y=max(0,abs(dy)-1)
                  do dx=-k_search,k_search
                     if (dz==0) then
                        if (dy<0 .or. (dy==0 .and. dx<=0)) cycle
                     endif
                     nx=cx1+dx
                     if (nx<1 .or. nx>nc(1)) cycle
                     gap_x=max(0,abs(dx)-1)
                     c2=nx+(ny-1)*nc(1)+(nz-1)*nc(1)*nc(2)
                     if (bound_two_level_pps(s)%cell_start(c2)>bound_two_level_pps(s)%cell_start(c2+1)-1) cycle
                     dmin2=cw2*dble(gap_x*gap_x+gap_y*gap_y+gap_z*gap_z)
                     rcut=max(task_max_a,bound_two_level_pps(s)%cell_max_a(c2))
                     if (dmin2>=rcut*rcut) cycle
                     do p1=p1_beg,p1_end
                        i1=bound_two_level_pps(s)%cell_particles(p1)
                        if (.not.keep(i1)) cycle
                        a1=a_table(max(1,min(4,ratios(i1))))
                        do p2=bound_two_level_pps(s)%cell_start(c2),bound_two_level_pps(s)%cell_start(c2+1)-1
                           i2=bound_two_level_pps(s)%cell_particles(p2)
                           if (.not.keep(i2)) cycle
                           a2=a_table(max(1,min(4,ratios(i2))))
                           dx_r=dble(q(1,i2))-dble(q(1,i1))
                           dy_r=dble(q(2,i2))-dble(q(2,i1))
                           dz_r=dble(q(3,i2))-dble(q(3,i1))
                           r2=dx_r*dx_r+dy_r*dy_r+dz_r*dz_r
                           r2_max_12=max(a1,a2)*max(a1,a2)
                           if (r2<r2_max_12) then
                              d=sqrt(r2)
                              phi_eps=bound_depth_continuous(d,eps)
                              if (r2<a1*a1) then
                                 delta_depth1=phi_eps-bound_depth_continuous(d,a1)
                                 dphi1=-real(factor*delta_depth1,4)
                                 bound_two_level_pps(s)%phi_team(i1,tid)=bound_two_level_pps(s)%phi_team(i1,tid)+dphi1
                              endif
                              if (r2<a2*a2) then
                                 delta_depth2=phi_eps-bound_depth_continuous(d,a2)
                                 dphi2=-real(factor*delta_depth2,4)
                                 bound_two_level_pps(s)%phi_team(i2,tid)=bound_two_level_pps(s)%phi_team(i2,tid)+dphi2
                              endif
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         !$omp end do
         !$omp do schedule(static)
         do i=1,np
            if (.not.keep(i)) cycle
            phi_sum=0.0
            do tid_k=1,nth
               phi_sum=phi_sum+bound_two_level_pps(s)%phi_team(i,tid_k)
            enddo
            phi_c2(i)=phi_sum
         enddo
         !$omp end do
         !$omp end parallel
      else
         do itask=1,ntasks
            c1=bound_two_level_pps(s)%task_c1(itask)
            p1_beg=bound_two_level_pps(s)%task_p1(itask)
            p1_end=bound_two_level_pps(s)%task_p2(itask)
            task_max_a=bound_two_level_pps(s)%task_max_a(itask)
            cz1=(c1-1)/(nc(1)*nc(2))+1
            cy1=mod(c1-1,nc(1)*nc(2))/nc(1)+1
            cx1=mod(mod(c1-1,nc(1)*nc(2)),nc(1))+1
            do p1=p1_beg,p1_end
               i1=bound_two_level_pps(s)%cell_particles(p1)
               if (.not.keep(i1)) cycle
               a1=a_table(max(1,min(4,ratios(i1))))
               do p2=p1+1,bound_two_level_pps(s)%cell_start(c1+1)-1
                  i2=bound_two_level_pps(s)%cell_particles(p2)
                  if (.not.keep(i2)) cycle
                  a2=a_table(max(1,min(4,ratios(i2))))
                  dx_r=dble(q(1,i2))-dble(q(1,i1))
                  dy_r=dble(q(2,i2))-dble(q(2,i1))
                  dz_r=dble(q(3,i2))-dble(q(3,i1))
                  r2=dx_r*dx_r+dy_r*dy_r+dz_r*dz_r
                  r2_max_12=max(a1,a2)*max(a1,a2)
                  if (r2<r2_max_12) then
                     d=sqrt(r2)
                     phi_eps=bound_depth_continuous(d,eps)
                     if (r2<a1*a1) then
                        delta_depth1=phi_eps-bound_depth_continuous(d,a1)
                        dphi1=-real(factor*delta_depth1,4)
                        phi_c2(i1)=phi_c2(i1)+dphi1
                     endif
                     if (r2<a2*a2) then
                        delta_depth2=phi_eps-bound_depth_continuous(d,a2)
                        dphi2=-real(factor*delta_depth2,4)
                        phi_c2(i2)=phi_c2(i2)+dphi2
                     endif
                  endif
               enddo
            enddo
            do dz=0,k_search
               nz=cz1+dz
               if (nz<1 .or. nz>nc(3)) cycle
               gap_z=max(0,abs(dz)-1)
               do dy=-k_search,k_search
                  ny=cy1+dy
                  if (ny<1 .or. ny>nc(2)) cycle
                  gap_y=max(0,abs(dy)-1)
                  do dx=-k_search,k_search
                     if (dz==0) then
                        if (dy<0 .or. (dy==0 .and. dx<=0)) cycle
                     endif
                     nx=cx1+dx
                     if (nx<1 .or. nx>nc(1)) cycle
                     gap_x=max(0,abs(dx)-1)
                     c2=nx+(ny-1)*nc(1)+(nz-1)*nc(1)*nc(2)
                     if (bound_two_level_pps(s)%cell_start(c2)>bound_two_level_pps(s)%cell_start(c2+1)-1) cycle
                     dmin2=cw2*dble(gap_x*gap_x+gap_y*gap_y+gap_z*gap_z)
                     rcut=max(task_max_a,bound_two_level_pps(s)%cell_max_a(c2))
                     if (dmin2>=rcut*rcut) cycle
                     do p1=p1_beg,p1_end
                        i1=bound_two_level_pps(s)%cell_particles(p1)
                        if (.not.keep(i1)) cycle
                        a1=a_table(max(1,min(4,ratios(i1))))
                        do p2=bound_two_level_pps(s)%cell_start(c2),bound_two_level_pps(s)%cell_start(c2+1)-1
                           i2=bound_two_level_pps(s)%cell_particles(p2)
                           if (.not.keep(i2)) cycle
                           a2=a_table(max(1,min(4,ratios(i2))))
                           dx_r=dble(q(1,i2))-dble(q(1,i1))
                           dy_r=dble(q(2,i2))-dble(q(2,i1))
                           dz_r=dble(q(3,i2))-dble(q(3,i1))
                           r2=dx_r*dx_r+dy_r*dy_r+dz_r*dz_r
                           r2_max_12=max(a1,a2)*max(a1,a2)
                           if (r2<r2_max_12) then
                              d=sqrt(r2)
                              phi_eps=bound_depth_continuous(d,eps)
                              if (r2<a1*a1) then
                                 delta_depth1=phi_eps-bound_depth_continuous(d,a1)
                                 dphi1=-real(factor*delta_depth1,4)
                                 phi_c2(i1)=phi_c2(i1)+dphi1
                              endif
                              if (r2<a2*a2) then
                                 delta_depth2=phi_eps-bound_depth_continuous(d,a2)
                                 dphi2=-real(factor*delta_depth2,4)
                                 phi_c2(i2)=phi_c2(i2)+dphi2
                              endif
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
   endsubroutine bound_two_level_pp_full

   subroutine bound_two_level_pp_update(s,np,q,keep,just_del_list,n_just_del, &
      spacing,app_val,a_val,ratios,phi_c2,nth)
      integer,intent(in) :: s,np,n_just_del,nth,just_del_list(n_just_del),ratios(np)
      real(4),intent(in) :: q(3,np),spacing,app_val
      logical,intent(in) :: keep(np)
      real(8),intent(in) :: a_val
      real(4),intent(inout) :: phi_c2(np)
      integer :: n_k,k,j,i,p,c_j,c2,cx_j,cy_j,cz_j,dx,dy,dz,nx,ny,nz,gap_x,gap_y,gap_z
      integer :: nc(3),k_search,tid,tid_k
      real(8) :: cw,cw2,dmin2,rcut,eps,factor,ddx,ddy,ddz,r2,d,ai,delta_depth,a_table(4)
      real(4) :: delta_phi,phi_sum
      if (nth<1 .or. nth>4) error stop 'two_level_pp_update: invalid team size'
      if (s<1 .or. s>fofcore .or. np<=0) error stop 'two_level_pp_update: invalid slot or count'
      if (spacing<=0.0 .or. app_val<=0.0) error stop 'two_level_pp_update: invalid length'
      n_k=count(keep(1:np))
      if (n_just_del<=0 .or. n_k<=0) return
      a_table(1)=3.5d0*dble(spacing)/1d0
      a_table(2)=bound_tl_split*dble(spacing)/2d0
      a_table(3)=0d0
      a_table(4)=bound_tl_split*dble(spacing)/4d0
      eps=dble(app_val)
      factor=a_val/(6d0*pi)
      nc=bound_two_level_pps(s)%ncell
      cw=bound_two_level_pps(s)%cw
      cw2=cw*cw
      k_search=bound_two_level_pps(s)%max_search_k
      if (nth>1) then
         if (.not.allocated(bound_two_level_pps(s)%phi_team)) then
            allocate(bound_two_level_pps(s)%phi_team(max(np,1),nth))
         else if (size(bound_two_level_pps(s)%phi_team,1)<np .or. &
                  size(bound_two_level_pps(s)%phi_team,2)<nth) then
            deallocate(bound_two_level_pps(s)%phi_team)
            allocate(bound_two_level_pps(s)%phi_team(max(np,1),nth))
         endif
         bound_two_level_pps(s)%phi_team(1:np,1:nth)=0.0
         !$omp parallel num_threads(nth) default(shared) &
         !$omp private(tid,tid_k,k,j,c_j,cz_j,cy_j,cx_j,dz,dy,dx,nz,ny,nx, &
         !$omp gap_z,gap_y,gap_x,c2,dmin2,rcut,p,i,ddx,ddy,ddz,r2,d,ai,delta_depth,delta_phi,phi_sum)
         !$omp single
         if (omp_get_num_threads()/=nth) error stop 'two_level_pp_update: actual team size mismatch'
         !$omp end single
         tid=omp_get_thread_num()+1
         !$omp do schedule(dynamic,1)
         do k=1,n_just_del
            j=just_del_list(k)
            c_j=bound_two_level_pps(s)%cell_id(j)
            cz_j=(c_j-1)/(nc(1)*nc(2))+1
            cy_j=mod(c_j-1,nc(1)*nc(2))/nc(1)+1
            cx_j=mod(mod(c_j-1,nc(1)*nc(2)),nc(1))+1
            do dz=-k_search,k_search
               nz=cz_j+dz
               if (nz<1 .or. nz>nc(3)) cycle
               gap_z=max(0,abs(dz)-1)
               do dy=-k_search,k_search
                  ny=cy_j+dy
                  if (ny<1 .or. ny>nc(2)) cycle
                  gap_y=max(0,abs(dy)-1)
                  do dx=-k_search,k_search
                     nx=cx_j+dx
                     if (nx<1 .or. nx>nc(1)) cycle
                     gap_x=max(0,abs(dx)-1)
                     c2=nx+(ny-1)*nc(1)+(nz-1)*nc(1)*nc(2)
                     if (bound_two_level_pps(s)%cell_start(c2)>bound_two_level_pps(s)%cell_start(c2+1)-1) cycle
                     dmin2=cw2*dble(gap_x*gap_x+gap_y*gap_y+gap_z*gap_z)
                     rcut=bound_two_level_pps(s)%cell_max_a(c2)
                     if (dmin2>=rcut*rcut) cycle
                     do p=bound_two_level_pps(s)%cell_start(c2),bound_two_level_pps(s)%cell_start(c2+1)-1
                        i=bound_two_level_pps(s)%cell_particles(p)
                        if (.not.keep(i)) cycle
                        ai=a_table(max(1,min(4,ratios(i))))
                        ddx=dble(q(1,i))-dble(q(1,j))
                        ddy=dble(q(2,i))-dble(q(2,j))
                        ddz=dble(q(3,i))-dble(q(3,j))
                        r2=ddx*ddx+ddy*ddy+ddz*ddz
                        if (r2<ai*ai) then
                           d=sqrt(r2)
                           delta_depth=bound_depth_continuous(d,eps)-bound_depth_continuous(d,ai)
                           delta_phi=-real(factor*delta_depth,4)
                           bound_two_level_pps(s)%phi_team(i,tid)=bound_two_level_pps(s)%phi_team(i,tid)-delta_phi
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
         !$omp end do
         !$omp do schedule(static)
         do i=1,np
            if (.not.keep(i)) cycle
            phi_sum=0.0
            do tid_k=1,nth
               phi_sum=phi_sum+bound_two_level_pps(s)%phi_team(i,tid_k)
            enddo
            if (phi_sum/=0.0) phi_c2(i)=phi_c2(i)+phi_sum
         enddo
         !$omp end do
         !$omp end parallel
      else
         do k=1,n_just_del
            j=just_del_list(k)
            c_j=bound_two_level_pps(s)%cell_id(j)
            cz_j=(c_j-1)/(nc(1)*nc(2))+1
            cy_j=mod(c_j-1,nc(1)*nc(2))/nc(1)+1
            cx_j=mod(mod(c_j-1,nc(1)*nc(2)),nc(1))+1
            do dz=-k_search,k_search
               nz=cz_j+dz
               if (nz<1 .or. nz>nc(3)) cycle
               gap_z=max(0,abs(dz)-1)
               do dy=-k_search,k_search
                  ny=cy_j+dy
                  if (ny<1 .or. ny>nc(2)) cycle
                  gap_y=max(0,abs(dy)-1)
                  do dx=-k_search,k_search
                     nx=cx_j+dx
                     if (nx<1 .or. nx>nc(1)) cycle
                     gap_x=max(0,abs(dx)-1)
                     c2=nx+(ny-1)*nc(1)+(nz-1)*nc(1)*nc(2)
                     if (bound_two_level_pps(s)%cell_start(c2)>bound_two_level_pps(s)%cell_start(c2+1)-1) cycle
                     dmin2=cw2*dble(gap_x*gap_x+gap_y*gap_y+gap_z*gap_z)
                     rcut=bound_two_level_pps(s)%cell_max_a(c2)
                     if (dmin2>=rcut*rcut) cycle
                     do p=bound_two_level_pps(s)%cell_start(c2),bound_two_level_pps(s)%cell_start(c2+1)-1
                        i=bound_two_level_pps(s)%cell_particles(p)
                        if (.not.keep(i)) cycle
                        ai=a_table(max(1,min(4,ratios(i))))
                        ddx=dble(q(1,i))-dble(q(1,j))
                        ddy=dble(q(2,i))-dble(q(2,j))
                        ddz=dble(q(3,i))-dble(q(3,j))
                        r2=ddx*ddx+ddy*ddy+ddz*ddz
                        if (r2<ai*ai) then
                           d=sqrt(r2)
                           delta_depth=bound_depth_continuous(d,eps)-bound_depth_continuous(d,ai)
                           delta_phi=-real(factor*delta_depth,4)
                           phi_c2(i)=phi_c2(i)-delta_phi
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
   endsubroutine bound_two_level_pp_update

   subroutine bound_two_level_pp_release
      integer :: s
      do s=1,fofcore
         if (allocated(bound_two_level_pps(s)%cell_id)) deallocate(bound_two_level_pps(s)%cell_id)
         if (allocated(bound_two_level_pps(s)%cell_particles)) deallocate(bound_two_level_pps(s)%cell_particles)
         if (allocated(bound_two_level_pps(s)%cell_start)) deallocate(bound_two_level_pps(s)%cell_start)
         if (allocated(bound_two_level_pps(s)%cell_max_a)) deallocate(bound_two_level_pps(s)%cell_max_a)
         if (allocated(bound_two_level_pps(s)%task_c1)) deallocate(bound_two_level_pps(s)%task_c1)
         if (allocated(bound_two_level_pps(s)%task_p1)) deallocate(bound_two_level_pps(s)%task_p1)
         if (allocated(bound_two_level_pps(s)%task_p2)) deallocate(bound_two_level_pps(s)%task_p2)
         if (allocated(bound_two_level_pps(s)%task_max_a)) deallocate(bound_two_level_pps(s)%task_max_a)
         if (allocated(bound_two_level_pps(s)%phi_team)) deallocate(bound_two_level_pps(s)%phi_team)
         bound_two_level_pps(s)%ncell=0; bound_two_level_pps(s)%ncell_total=0
         bound_two_level_pps(s)%u_min=0d0; bound_two_level_pps(s)%cw=0d0
         bound_two_level_pps(s)%max_search_k=1; bound_two_level_pps(s)%max_cutoff=0d0
         bound_two_level_pps(s)%ntasks=0
      enddo
   endsubroutine bound_two_level_pp_release

   ! FoF契约：每个CRS区间恰有np_list个互异合法内部下标；稳定ID唯一。
   ! N 已满足 np_halo_min，xv 为有限的六维相空间且坐标已展开周期边界。
   ! 调用期间 xv/CRS/pid 不可修改；DSU构建保证内部成员区间合法且互不重叠。
   ! 外层设置 OMP_DYNAMIC=FALSE、OMP_PLACES=cores、OMP_PROC_BIND=close；
   ! 这只是运行约定，模块不验证物理 CPU affinity。仅链接单线程 FFTW3f。
   ! 批入口必须由外层串行调用，不可同时或嵌套调用；当前共享调度/资源不支持重入。
   ! Config 6/8/9 selects a 16/32/64 CPU budget; default is 32 physical cores.
   ! final_members 为候选容量前缀下的原xv内部下标；PID只用于稳定排序。
   subroutine bound_check_batch(head,np_list,x_mean,v_mean,np_bound,x_bound,v_bound,status,configuration, &
       final_members)
      integer,intent(in) :: head(:),np_list(:)
      real(4),intent(in) :: x_mean(:,:),v_mean(:,:)
      integer,intent(out) :: np_bound(:),status(:)
      real(4),intent(out) :: x_bound(:,:),v_bound(:,:)
      integer,optional,intent(in) :: configuration
      integer,optional,intent(out) :: final_members(:)
      integer,allocatable :: result_offset(:)
      integer(8) :: result_size
      type(bound_schedule_type) :: tasks
      integer :: config,g,l,s,n,cap,pcap,ios,st,h,pos,next(4),iphase,nslots,budget
      integer :: phase_workers(4,2),phase_threads(2),i,cnt,nphases
      logical,allocatable :: halo_processed(:)
      real(8) :: phase_start,batch_start

       batch_start=omp_get_wtime()
       config=8
       if (present(configuration)) config=configuration
      ! 当前每个新批先释放旧工作槽，因此 FFT 计划按批重建；各类只读 Green 核仍保留复用。
       call bound_release
       phase_workers=0; phase_threads=1; nphases=2
       select case(config)
        case(6)
          budget=16
          phase_workers(:,1)=[0,0,0,8]
          phase_workers(:,2)=[3,5,8,0]
        case(8)
          budget=32
          phase_workers(:,1)=[0,0,0,16]
          phase_workers(:,2)=[6,10,16,0]
        case(9)
          budget=64
          phase_workers(:,1)=[0,0,0,32]
          phase_workers(:,2)=[12,20,32,0]
        case default
          error stop 'bound_check_batch: configuration must be 6, 8 or 9'
        end select
        if (nphases==2) phase_threads=[2,1]
       if (ngrid /= 4) error stop 'bound_check_batch: fixed groups require ng >= 256'
       if (omp_in_parallel()) error stop 'bound_check_batch: caller must be serial'
       if (size(np_list)/=size(head) .or. size(np_bound)/=size(head) .or. size(status)/=size(head)) &
          error stop 'bound_check_batch: incompatible output sizes'
       if (size(x_mean,1)/=3 .or. size(v_mean,1)/=3 .or. size(x_bound,1)/=3 .or. size(v_bound,1)/=3) &
          error stop 'bound_check_batch: expected three-vector arrays'
       if (size(x_mean,2)/=size(head) .or. size(v_mean,2)/=size(head) .or. &
           size(x_bound,2)/=size(head) .or. size(v_bound,2)/=size(head)) &
          error stop 'bound_check_batch: incompatible vector counts'
       if (.not.ieee_is_finite(a_bound) .or. a_bound<=0d0 .or. .not.ieee_is_finite(b_bound)) &
          error stop 'bound_check_batch: invalid physical context'
       if (present(final_members)) then
          if (.not.allocated(bound_member_table)) error stop 'bound_check_batch: final members requires CRS input'
          result_size=sum(int(np_list,8))
          if (result_size>int(huge(0),8)) error stop 'bound_check_batch: compact member output too large'
          if (size(final_members,kind=8)/=result_size) error stop 'bound_check_batch: final members size mismatch'
          allocate(result_offset(size(head)))
          result_size=1_8
          do h=1,size(head)
             result_offset(h)=int(result_size)
             result_size=result_size+int(np_list(h),8)
          enddo
          final_members=0
       endif
       call bound_classify(head,np_list,budget,int(ng),tasks)
       allocate(halo_processed(size(head)),stat=ios)
       if (ios/=0) error stop 'bound_check_batch: processed allocation failed'
        halo_processed=.false.
       next=tasks%first
       np_bound=0; x_bound=0.0; v_bound=0.0; status=-1
        call omp_set_max_active_levels(2)
        do iphase=1,nphases
       phase_start=omp_get_wtime()
       bound_workers=phase_workers(:,iphase)
       bound_c2_threads=phase_threads(iphase)
       nslots=sum(bound_workers)
       bound_tl_slot_count=max(1,bound_workers(4)*bound_c2_threads)
       if (nslots>fofcore .or. bound_tl_slot_count>bound_tl_max_slots) &
          error stop 'bound_check_batch: slot budget exceeded'
       write(*,'(A,I0,A,I0,A,I0)') '[SCHEDULE] phase=',iphase,' owners=',nslots,' team_threads=',bound_c2_threads
       worker_group=0; worker_slot=0
       s=0
      do g=4,1,-1
         do l=1,bound_workers(g)
            s=s+1
            worker_group(s)=g; worker_slot(s)=l
         enddo
      enddo
       allocate(bound_mesh(nslots),bound_mem(nslots),stat=ios)
      if (ios/=0) error stop 'bound_check_batch: worker allocation failed'
       do g=1,4
          if (tasks%pm_capacity(g)==0 .or. bound_workers(g)==0) cycle
         n=bn(g)
         allocate(bound_storage(g)%rho(n+2,n,n,bound_workers(g)),stat=ios)
         if (ios/=0) error stop 'bound_check_batch: PM allocation failed'
         ! n 为偶数：实数块 (n+2)*n*n 与复数块 (n/2+1)*n*n 字节数相同。
         call c_f_pointer(c_loc(bound_storage(g)%rho),bound_storage(g)%ft, &
            [n/2+1,n,n,bound_workers(g)])
      enddo
       do s=1,nslots
         g=worker_group(s); cap=tasks%capacity(g); pcap=tasks%pm_capacity(g)
         if (cap==0) cycle
         call bound_init_rho(s,merge(g,0,pcap>0),ng,cap,st)
         if (st/=0) error stop 'bound_check_batch: slot initialization failed'
         ! 每槽按本组最大 N 预留公共缓冲；PM 专用缓冲按本组最大 PM 候选数 pcap 分配。
         associate(b=>bound_mem(s))
            allocate(b%order(cap),b%q(3,cap),b%phi(cap),b%energy(cap), &
               b%just_del_list(cap),b%phi_direct(cap),stat=ios)
            if (ios/=0) error stop 'bound_check_batch: particle allocation failed'
             allocate(b%self(pcap),b%first_tsc(3,pcap),b%wt_tsc(0:2,3,pcap), &
                b%phi_c1(pcap),b%phi_c2(pcap),stat=ios)
             if (ios/=0) error stop 'bound_check_batch: PM particle allocation failed'
             allocate(b%phi_c2_team(pcap,bound_c2_max_threads),stat=ios)
             if (ios/=0) error stop 'bound_check_batch: team buffer allocation failed'
         endassociate
         if (pcap>0) then
            allocate(bound_cells(s)%cell_id(pcap),bound_cells(s)%cell_particles(pcap),stat=ios)
            if (ios/=0) error stop 'bound_check_batch: cell particle allocation failed'
            call bound_init_gk(s,st)
         endif
       enddo
       if (bound_workers(4)>0) call bound_two_level_prepare_schedule(tasks)
      ! next(g) 是第 g 组的下一个队列位置；空组 first>last，所属 worker 会直接退出。
        !$omp parallel num_threads(nslots) default(shared) &
        !$omp private(s,g,pos,h,i,cnt)
      !$omp single
        if (omp_get_num_threads()/=nslots) error stop 'bound_check_batch: actual team size mismatch'
      !$omp end single
      s=omp_get_thread_num()+1
       g=worker_group(s)
       do
         ! 原子领取的是队列位置，不是粒子索引；同组各 worker 不会取得同一 Halo。
          !$omp atomic capture
          pos=next(g)
          next(g)=next(g)+1
          !$omp end atomic
          if (pos>tasks%last(g)) exit
          h=tasks%queue(pos)
         ! 只装载已选好的任务参数，求解器不再重选网格；直接任务方法=0，但仍归48工作组。
         bound_mem(s)%method_gid=tasks%method_gid(h)
         bound_mem(s)%anchor=tasks%anchor(:,h)
         if (tasks%method_gid(h) >= 1) then
            bound_mem(s)%origin=tasks%origin(:,h)
         else
            bound_mem(s)%origin=0.0
         endif
         bound_mesh(s)%level=tasks%level(h)
         bound_mesh(s)%spacing=tasks%spacing(h)
         bound_mesh(s)%span=dble(tasks%span(h))
         ! h 始终为原输入任务编号；各线程写不同的 h，输出顺序不随领取顺序改变。
          call bound_check(h,np_list(h),x_mean(:,h),v_mean(:,h), &
              np_bound(h),x_bound(:,h),v_bound(:,h),status(h),slot=s)
          if (present(final_members)) then
             if (status(h)==0 .or. status(h)==5) then
                cnt=0
                do i=1,np_list(h)
                   if (.not.bound_mem(s)%keep(i)) cycle
                   cnt=cnt+1
                   final_members(result_offset(h)+cnt-1)=bound_mem(s)%member(i)
                enddo
                if (cnt/=np_bound(h)) error stop 'bound_check_batch: retained member count mismatch'
             endif
          endif
          halo_processed(h)=.true.
       enddo
       !$omp end parallel
       call bound_release
       write(*,'(A,I0,A,F12.4)') '[BOUND_PHASE] phase=',iphase,' total_sec=',omp_get_wtime()-phase_start
       enddo
       if (.not.all(halo_processed)) error stop 'bound_check_batch: unprocessed candidates'
       if (any(status/=0 .and. status/=1 .and. status/=5)) error stop 'bound_check_batch: unexpected status'
       write(*,'(A,3(I0,1X),A,F12.4)') '[BOUND_STATUS] OK/TOO_FEW/MAX_ITER=', &
          count(status==0),count(status==1),count(status==5),' total_sec=',omp_get_wtime()-batch_start
    endsubroutine bound_check_batch

   ! 串行销毁所有槽的 FFT 计划，再释放原位工作区和成员/CSR缓存；只读 Green 核不在此销毁。
   subroutine bound_release
      integer :: s,g,st
      if (allocated(bound_mesh)) then
         do s=1,size(bound_mesh)
            call bound_init_rho(s,0,ng,0,st)
         enddo
         deallocate(bound_mesh,bound_mem)
      endif
      do g=1,4
         nullify(bound_storage(g)%ft)
         if (associated(bound_storage(g)%rho)) deallocate(bound_storage(g)%rho)
      enddo
      do s=1,fofcore
         if (allocated(bound_cells(s)%cell_id)) deallocate(bound_cells(s)%cell_id)
         if (allocated(bound_cells(s)%cell_start)) deallocate(bound_cells(s)%cell_start)
          if (allocated(bound_cells(s)%cell_particles)) deallocate(bound_cells(s)%cell_particles)
          if (allocated(bound_cells(s)%task_c1)) deallocate(bound_cells(s)%task_c1)
          if (allocated(bound_cells(s)%task_p1)) deallocate(bound_cells(s)%task_p1)
          if (allocated(bound_cells(s)%task_p2)) deallocate(bound_cells(s)%task_p2)
          bound_cells(s)%ncell=0; bound_cells(s)%ncell_total=0; bound_cells(s)%ntasks=0
          call bound_two_level_release_halo(s)
       enddo
       call bound_release_inner_slots
       call bound_two_level_pp_release
   endsubroutine bound_release

   !----------------------------------------------------------------------
   ! bound_grid_levels: 返回指定网格类别启用的 spacing 档位表（按覆盖升序）
   !
   ! 规则（与 TEST_PLAN.md 动态多分辨率网格一致）：
   !   1. spacing 取自几何序列 1/32, 1/16, 1/8, 1/4, 1/2, ...；
   !   2. 只保留物理覆盖 m*spacing 严格小于 15 Mpc 对应跨度 S_15 的档；
   !   3. 最大启用网格（gid == ngrid）再追加精确 spacing_max，使其有效覆盖
   !      严格等于 15 Mpc，不加物理覆盖余量；
   !   4. 返回表已按物理覆盖升序排列，相邻档覆盖连续，不出现分辨率断崖。
   ! 未启用（n_grid > ng）的网格类别返回 nlev = 0。
   ! sp_out 为调用方提供的缓冲区，长度不得小于实际档数（上限 nlev_max）。
   !----------------------------------------------------------------------
   subroutine bound_grid_levels(gid, nlev, sp_out)
      implicit none
      integer,intent(in) :: gid
      integer,intent(out) :: nlev
      real(8),intent(out) :: sp_out(:)
      integer m
      real(8) sp

      nlev = 0
      if (gid < 1 .or. gid > ngrid) return
      if (size(sp_out) < 1) return
      m = bn(gid)/2 - 3
      sp = 1d0/32d0
      do while (dble(m)*sp < s_15 .and. nlev < size(sp_out))
         nlev = nlev + 1
         sp_out(nlev) = sp
         sp = 2d0*sp
      enddo
      if (gid == ngrid .and. nlev < size(sp_out)) then
         nlev = nlev + 1
         sp_out(nlev) = spacing_max
      endif
   endsubroutine bound_grid_levels

   !----------------------------------------------------------------------
   ! bound_build_ladder: 惰性建立全局档位阶梯（各网格全部启用档按覆盖升序）
   ! 只在首次调用时构建，之后重用；档数由启用网格与 spacing 上限决定，使用 nlev_max 缓冲。
   !----------------------------------------------------------------------
   subroutine bound_build_ladder
      implicit none
      integer g,l,i,j,k,nlev,nt
      real(8) sp(32),ts(nlev_max),tc(nlev_max),tmp
      integer tg(nlev_max),tl(nlev_max),itmp

      if (bound_ntotal > 0) return
      nt = 0
      do g = 1, ngrid
         call bound_grid_levels(g, nlev, sp)
         do l = 1, nlev
            nt = nt + 1
            tg(nt) = g; tl(nt) = l; ts(nt) = sp(l)
            tc(nt) = dble(bn(g)/2 - 3)*sp(l)
         enddo
      enddo
      ! 按物理覆盖升序做选择排序（档数很少，串行执行一次）
      do j = 1, nt - 1
         k = j
         do i = j + 1, nt
            if (tc(i) < tc(k)) k = i
         enddo
         if (k /= j) then
            itmp = tg(j); tg(j) = tg(k); tg(k) = itmp
            itmp = tl(j); tl(j) = tl(k); tl(k) = itmp
            tmp = tc(j); tc(j) = tc(k); tc(k) = tmp
            tmp = ts(j); ts(j) = ts(k); ts(k) = tmp
         endif
      enddo
      bound_tab_gid(1:nt) = tg(1:nt)
      bound_tab_lev(1:nt) = tl(1:nt)
      bound_tab_sp(1:nt) = ts(1:nt)
      bound_tab_cov(1:nt) = tc(1:nt)
      bound_ntotal = nt
   endsubroutine bound_build_ladder

   !----------------------------------------------------------------------
   ! bound_fgl_init: 惰性初始化固定选档候选序列与阶梯镜像
   !----------------------------------------------------------------------
   subroutine bound_fgl_init
      implicit none
      integer :: gid, nlev
      real(8) :: sp(bound_fgl_maxlev)

      if (bound_fgl_ready) return
      do gid = 1, ngrid
         call bound_grid_levels(gid, nlev, sp)
         bound_fgl_nlev(gid) = nlev
         if (nlev > 0) then
            bound_fgl_sp(gid, 1:nlev) = sp(1:nlev)
            bound_fgl_cov(gid, 1:nlev) = dble(bn(gid)/2 - 3) * sp(1:nlev)
         endif
      enddo

      ! band 1 (1000 < N <= 5000): 48/0.125, 64/0.125, 64/0.25, 128(0.125,0.25,0.5), 256
      bound_fgl_ncand(1) = 11
      bound_fgl_cgid(1, 1:11) = [1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4]
      bound_fgl_csp(1, 1:11) = [0.125d0, 0.125d0, 0.25d0, 0.125d0, 0.25d0, 0.5d0, &
         0.0625d0, 0.125d0, 0.25d0, 0.5d0, bound_fgl_term_sentinel]

      ! band 2 (5000 < N <= 20000): 128/0.0625, 128/0.125, D>7.625 直接进入 256 序列
      bound_fgl_ncand(2) = 10
      bound_fgl_cgid(2, 1:10) = [3, 3, 4, 4, 4, 4, 4, 4, 4, 4]
      bound_fgl_csp(2, 1:10) = [0.0625d0, 0.125d0, 0.0625d0, 0.125d0, 0.25d0, 0.5d0, &
         1.0d0, 2.0d0, 4.0d0, bound_fgl_term_sentinel]

      ! band 3 (N > 20000): 256 跨度表
      bound_fgl_ncand(3) = 8
      bound_fgl_cgid(3, 1:8) = [4, 4, 4, 4, 4, 4, 4, 4]
      bound_fgl_csp(3, 1:8) = [0.0625d0, 0.125d0, 0.25d0, 0.5d0, &
         1.0d0, 2.0d0, 4.0d0, bound_fgl_term_sentinel]

      bound_fgl_ready = .true.
   endsubroutine bound_fgl_init

   !----------------------------------------------------------------------
   ! bound_select_grid_and_level: 正式统一选档子程序
   ! 输入：nParticles（初始候选粒子数）、span（展开后真实最大轴向跨度）
   ! 输出：gid_out（0=直接求和, 1..4=PM网格）、lev_out、spacing_out、status
   ! 规则来源：FIXED_GRID_CORE_ALLOCATION.md
   !----------------------------------------------------------------------
   subroutine bound_select_grid_and_level(nParticles, span, gid_out, lev_out, spacing_out, status)
      implicit none
      integer(8),intent(in) :: nParticles
      real(8),intent(in) :: span
      integer,intent(out) :: gid_out, lev_out, status
      real(8),intent(out) :: spacing_out
      integer :: band, slot, g, l, lev
      real(8) :: sp

      call bound_fgl_init
      gid_out = 0; lev_out = 0; spacing_out = 0d0; status = BOUND_SEL_OK

      ! 初始 N 已达到 np_halo_min，由提交任务的 FoF 入口保证。
      if (nParticles <= 1000_8) then
         gid_out = BOUND_METHOD_DIRECT
         status = BOUND_SEL_OK
         return
      endif
      if (nParticles <= 5000_8) then
         band = 1
      else if (nParticles <= 20000_8) then
         band = 2
      else
         band = 3
      endif

      do slot = 1, bound_fgl_ncand(band)
         g = bound_fgl_cgid(band, slot)
         if (g < 1 .or. g > ngrid) cycle
         sp = bound_fgl_csp(band, slot)
         if (sp == bound_fgl_term_sentinel) sp = spacing_max
         lev = 0
         do l = 1, bound_fgl_nlev(g)
            if (abs(bound_fgl_sp(g, l) - sp) <= 1d-12 * max(1d0, abs(sp))) then
               lev = l
               exit
            endif
         enddo
         if (lev == 0) cycle
         if (bound_fgl_cov(g, lev) < span) cycle
         gid_out = g
         lev_out = lev
         spacing_out = bound_fgl_sp(g, lev)
         return
      enddo
      ! 候选序列没有可用档时，查询当前启用阶梯的覆盖档；固定批入口要求 ng>=256。
      if (span <= s_15) then
         call bound_select_level(span, g, lev)
         if (g >= 1 .and. g <= ngrid .and. lev >= 1) then
            gid_out = g
            lev_out = lev
            spacing_out = bound_fgl_sp(g, lev)
            return
         endif
      endif
      status = BOUND_SEL_NO_COMBO
   endsubroutine bound_select_grid_and_level

   !----------------------------------------------------------------------
   ! bound_select_level: 按全局档位阶梯返回覆盖该跨度的最小覆盖 PM 档
   ! （仅供阶梯覆盖查询/诊断对照，正式自动选档统一使用 bound_select_grid_and_level）
   !----------------------------------------------------------------------
   subroutine bound_select_level(span, gsel, lsel)
      implicit none
      real(8),intent(in) :: span
      integer,intent(out) :: gsel,lsel
      integer i

      gsel = 0; lsel = 0
      call bound_build_ladder
      do i = 1, bound_ntotal
         if (bound_tab_cov(i) >= span) then
            gsel = bound_tab_gid(i)
            lsel = bound_tab_lev(i)
            return
         endif
      enddo
   endsubroutine bound_select_level

   !----------------------------------------------------------------------
   ! bound_grid_reject: 网格选择阶段的统一终止出口（不新增返回状态码）
   !   · 候选真实展开跨度大于 15 Mpc（D_max）时使用；
   !   · 或调用方选定的工作槽网格无法容纳该跨度（档位与调度不一致）时使用。
   ! 打印候选标识、实际跨度（标准格点与共动 Mpc 双单位）、15 Mpc 上限、S_15、
   ! box/h0/ng/nn/ng_global、n_max 与 spacing_max，随后以非零状态终止。
   ! 不增大网格、不放粗 spacing、不截断成员、不写任何有效结果。
   !----------------------------------------------------------------------
   subroutine bound_grid_reject(where, ident, span, reason)
      implicit none
      character(*),intent(in) :: where,reason
      integer,intent(in) :: ident
      real(8),intent(in) :: span
      real(8) l_ng

      l_ng = dble(box)/(dble(h0)*dble(ng_global))
      print*,'==============================================================='
      print*,'bound grid selection: candidate rejected, no valid result written'
      print*,'  subroutine           = ',trim(where)
      print*,'  reason               = ',trim(reason)
      print*,'  candidate identifier = ',ident
      print*,'  actual span          = ',span,' standard grid cells'
      print*,'  actual span          = ',span*l_ng,' Mpc (comoving)'
      print*,'  upper limit D_max    = 15 Mpc = ',s_15,' standard grid cells'
      print*,'  box [Mpc/h], h0      = ',box,h0
      print*,'  ng, nn, ng_global    = ',ng,nn,ng_global
      print*,'  n_max, m_max         = ',n_max,m_max
      print*,'  spacing_max          = ',spacing_max
      print*,'  nlev ladder size     = ',bound_ntotal,' ngrid = ',ngrid
      print*,'  L_ng [Mpc]           = ',l_ng
      print*,'==============================================================='
      error stop 1
   endsubroutine bound_grid_reject

   !----------------------------------------------------------------------
   ! bound_init_rho: 初始化单个 worker 的 FFT 计划、成员表与 keep 缓冲；或释放这些资源
   !
   ! 参数说明：
   !   s           [in]  : 当前工作槽编号 (1 <= s <= fofcore)
   !   grid_id     [in]  : 网格类别 (0=直接求和, 1=48, 2=64, 3=128, 4=256)
   !   ng_in       [in]  : 保留原接口；当前网格尺寸由模块参数 bn 确定
   !   np_capacity [in]  : 该槽预估最大成员粒子容量（np_capacity == 0 表示释放该槽资源）
   !   status      [out] : 返回状态码（0 表示成功，非 0 表示参数错误或分配失败）
   !
   ! 生命周期与并发约定：
   !   1. 必须在并行区外串行调用；
   !   2. 预分配固定内存空间，避免循环中频繁 allocate/deallocate 造成内存碎片；
   !   3. np_capacity == 0 时调用 sfftw_destroy_plan 销毁计划，并释放 member/keep；其余缓冲由批释放处理。
   !----------------------------------------------------------------------
   subroutine bound_init_rho(s, grid_id, ng_in, np_capacity, status)
      implicit none
      integer,intent(in) :: s,grid_id,np_capacity
      integer(8),intent(in) :: ng_in
      integer,intent(out) :: status
      integer n,ios,l

      status = 3
      if (s < 1 .or. s > fofcore) error stop 'bound_init_rho: invalid worker'
      if (.not.allocated(bound_mesh) .or. .not.allocated(bound_mem)) &
         error stop 'bound_init_rho: resources not initialized'
      l=worker_slot(s)

      !------------------------------------------------------------------
      ! 分支 1: np_capacity == 0 触发显式释放逻辑
      !------------------------------------------------------------------
      if (np_capacity == 0) then
         if (bound_mesh(s)%plan /= 0) then
            call sfftw_destroy_plan(bound_mesh(s)%plan)
            bound_mesh(s)%plan = 0
         endif
         if (bound_mesh(s)%iplan /= 0) then
            call sfftw_destroy_plan(bound_mesh(s)%iplan)
            bound_mesh(s)%iplan = 0
         endif
         if (allocated(bound_mem(s)%member)) deallocate(bound_mem(s)%member)
         if (allocated(bound_mem(s)%keep))   deallocate(bound_mem(s)%keep)
         bound_mesh(s)%grid_id = 0
         bound_mesh(s)%level = 0
         bound_mesh(s)%required_grid = 0
         bound_mesh(s)%spacing = 0d0
         bound_mesh(s)%span = 0d0
         status = 0
         return
      endif

      !------------------------------------------------------------------
      ! 分支 2: 初始化当前 worker 对应的组内工作区切片及 FFT 计划
      ! 只接受实际启用的网格类别：ng<256 时第 4 槽被禁用，不会创建计划。
      !------------------------------------------------------------------
      if (grid_id < 0 .or. grid_id > ngrid) error stop 'bound_init_rho: invalid grid'
      if (np_capacity < 1) error stop 'bound_init_rho: invalid capacity'

      ! 若槽内已有旧计划，先安全清理
      if (bound_mesh(s)%plan /= 0) then
         call sfftw_destroy_plan(bound_mesh(s)%plan)
         bound_mesh(s)%plan = 0
      endif
      if (bound_mesh(s)%iplan /= 0) then
         call sfftw_destroy_plan(bound_mesh(s)%iplan)
         bound_mesh(s)%iplan = 0
      endif

      ! 创建单精度三维实数-复数 FFTW 原位变换计划 (64 = FFTW_ESTIMATE 模式)
      if (grid_id >= 1) then
         n = bn(grid_id)
         call sfftw_plan_dft_r2c_3d(bound_mesh(s)%plan,n,n,n, &
            bound_storage(grid_id)%rho(:,:,:,l),bound_storage(grid_id)%ft(:,:,:,l),64)
         call sfftw_plan_dft_c2r_3d(bound_mesh(s)%iplan,n,n,n, &
            bound_storage(grid_id)%ft(:,:,:,l),bound_storage(grid_id)%rho(:,:,:,l),64)
         if (bound_mesh(s)%plan == 0 .or. bound_mesh(s)%iplan == 0) &
            error stop 'bound_init_rho: FFT plan creation failed'
      endif

      ! 分配独立的成员粒子索引与束缚状态掩码缓存
      if (allocated(bound_mem(s)%member)) deallocate(bound_mem(s)%member)
      if (allocated(bound_mem(s)%keep))   deallocate(bound_mem(s)%keep)
      allocate(bound_mem(s)%member(np_capacity), bound_mem(s)%keep(np_capacity), stat=ios)
      if (ios /= 0) error stop 'bound_init_rho: member allocation failed'

      bound_mesh(s)%grid_id = grid_id
      bound_mesh(s)%level = 0
      bound_mesh(s)%required_grid = 0
      bound_mesh(s)%spacing = 0d0
      bound_mesh(s)%span = 0d0
      status = 0
   endsubroutine bound_init_rho

   !----------------------------------------------------------------------
   ! bound_init_gk: 预计算指定网格类别的频域 Green 核与自身势基核
   !
   ! 参数说明：
   !   s      [in]  : 已完成 bound_init_rho 初始化的工作槽编号
   !   status [out] : 返回状态码（0 表示就绪，非 0 表示异常）
   !
   ! 物理与实现细节：
   !   1. C1 核在局部网格单位下生成，分裂软化固定为 h_split=3.5，而不是 app/spacing；
   !   2. 目标软化 app/spacing 在 C2 残差中使用，补回连续软化势与 C1 平滑势之差；
   !   3. green_phi 直接在频域生成；逆 FFT 后提取 5x5x5 离散响应，用于 TSC 自身势扣除；
   !   4. 同网格类别的 Green 核对所有线程只读共享，只在首次遇到该类别时计算一次。
   !----------------------------------------------------------------------
   subroutine bound_init_gk(s, status)
      implicit none
      integer,intent(in) :: s
      integer,intent(out) :: status
      integer gid

      status = 3
      if (s < 1 .or. s > size(bound_mesh)) error stop 'bound_init_gk: invalid worker'
      gid = bound_mesh(s)%grid_id
      if (gid == 0) then
         status = 0
         return
      endif
      if (gid < 1 .or. gid > ngrid) error stop 'bound_init_gk: invalid grid'

      call bound_init_gk_grid(s,gid,bound_storage(gid)%rho(:,:,:,worker_slot(s)), &
         bound_storage(gid)%ft(:,:,:,worker_slot(s)),status)
      if (status /= 0) error stop 'bound_init_gk: resource initialization failed'
   endsubroutine bound_init_gk

   subroutine bound_init_gk_grid(s, gid, rho, ft, status)
      implicit none
      integer,intent(in) :: s,gid
      real(4),intent(inout) :: rho(:,:,:)
      complex(4),intent(inout) :: ft(:,:,:)
      integer,intent(out) :: status
      integer n,nlev,ix,iy,iz,i,j,k,ia,ja,ka,ios
      real(8) :: sp_levels(32)
      real(8) :: kval(3), kalias(3), kmag, uk2, asum, bsum
      real(8) :: ssoft, srange, target_phi, cutoff_factor
      real(8) :: y, kr, radius
      real(4) :: apm_soft, apm_range
      integer :: assignment_order

      status = 3
      n = bn(gid)

      ! 按动态多分辨率网格规则取得该网格类别启用的 spacing 档位表。
      ! n_grid > ng 的网格类别被禁用（nlev = 0），不生成核也不创建计划。
      call bound_grid_levels(gid, nlev, sp_levels)
      if (nlev < 1) then
         print*,'bound_init_gk_grid: grid class not enabled (n_grid > ng), gid =',gid,' ng =',ng
         return
      endif

      ! 若当前网格类别的核已就绪且尺寸、档位数、软化、ng 全部相符，直接复用，避免重复计算
      if (allocated(bound_gk(gid)%green_phi) .or. allocated(bound_gk(gid)%spacing)) then
         if (allocated(bound_gk(gid)%green_phi) .and. allocated(bound_gk(gid)%spacing) .and. &
            bound_gk(gid)%n == n .and. bound_gk(gid)%nlevel == nlev .and. &
            abs(bound_gk(gid)%softening - dble(app)) < 1d-12 .and. &
            bound_gk(gid)%ng_standard == ng) then
            status = 0
            return
         else
            if (allocated(bound_gk(gid)%spacing)) deallocate(bound_gk(gid)%spacing)
            if (allocated(bound_gk(gid)%green_phi)) deallocate(bound_gk(gid)%green_phi)
         endif
      endif

      ! 1. 保存合法 spacing 等元数据
      bound_gk(gid)%n = n
      bound_gk(gid)%nlevel = nlev
      bound_gk(gid)%ng_standard = int(ng)
      bound_gk(gid)%softening = dble(app)
      allocate(bound_gk(gid)%spacing(nlev), stat=ios)
      if (ios /= 0) error stop 'bound_init_gk_grid: spacing allocation failed'
      bound_gk(gid)%spacing = sp_levels(1:nlev)

      ! 2. 生成当前 C1 最优频域势核 green_phi
      allocate(bound_gk(gid)%green_phi(n/2+1, n, n), stat=ios)
      if (ios /= 0) error stop 'bound_init_gk_grid: Green allocation failed'

       call bound_calc_green_and_self(n,bound_h_split,0d0,bound_mesh(s)%iplan, &
          rho,ft,bound_gk(gid)%green_phi,bound_gk(gid)%self_phi)

      status = 0
    endsubroutine bound_init_gk_grid

   ! Preserve the frozen alias sum and per-node arithmetic, with Green16 setup.
   subroutine bound_calc_green_and_self(n,apm_soft,apm_range,iplan,rho_buf,ft_buf,green_out,self_out,truncated)
      integer,intent(in) :: n
      real(8),intent(in) :: apm_soft,apm_range
      integer(8),intent(in) :: iplan
      real(4),intent(inout) :: rho_buf(:,:,:)
      complex(4),intent(inout) :: ft_buf(:,:,:)
      complex(4),intent(out) :: green_out(:,:,:)
      real(4),intent(out) :: self_out(-2:2,-2:2,-2:2)
      logical,optional,intent(in) :: truncated
      logical :: is_truncated
      integer :: ix,iy,iz,i,j,k,ia,ja,ka,assignment_order,eff_threads
      real(8) :: kval(3),kalias(3),kmag,uk2,asum,bsum,ssoft,srange,target_phi,cutoff_factor,y,kr,radius
      is_truncated=.true.
      if (present(truncated)) is_truncated=truncated
      assignment_order=2
      radius=0.5d0*dble(n)
      eff_threads=1
      if (.not.omp_in_parallel()) eff_threads=bound_green_threads
      !$omp parallel do collapse(2) schedule(static) num_threads(eff_threads) &
      !$omp private(k,j,i,kval,asum,bsum,ka,ja,ia,kalias,kmag,uk2,y,ssoft,srange,kr,cutoff_factor,target_phi)
      do k=1,n
         do j=1,n
            do i=1,n/2+1
               kval=[bound_fft_wave_number(i,n),bound_fft_wave_number(j,n),bound_fft_wave_number(k,n)]
               asum=0d0; bsum=0d0
               do ka=-n_int,n_int
                  do ja=-n_int,n_int
                     do ia=-n_int,n_int
                        kalias=kval+2d0*pi*dble([ia,ja,ka])
                        kmag=sqrt(sum(kalias*kalias))
                        uk2=bound_assignment_window_sq(kalias,assignment_order)
                        asum=asum+uk2
                        if (kmag==0d0) then
                           if (apm_range==0d0) then
                              target_phi=-4d0*pi*radius*radius/6d0
                           else
                              if (is_truncated) then
                                 target_phi=0d0
                              else
                                 target_phi=-2d0*pi*(apm_range**2-apm_soft**2)/15d0
                              endif
                           endif
                        else
                           y=apm_soft*kmag/2d0
                           ssoft=12d0*y**(-4)*(2d0-2d0*cos(y)-y*sin(y))
                           if (apm_range==0d0) then
                              srange=0d0
                           else
                              y=apm_range*kmag/2d0
                              srange=12d0*y**(-4)*(2d0-2d0*cos(y)-y*sin(y))
                           endif
                           if (is_truncated) then
                              kr=kmag*radius
                              cutoff_factor=1d0-sin(kr)/kr
                           else
                              cutoff_factor=1d0
                           endif
                           target_phi=-4d0*pi*(ssoft*ssoft-srange*srange)*cutoff_factor/(kmag*kmag)
                        endif
                        bsum=bsum+uk2*target_phi
                     enddo
                  enddo
               enddo
               if (asum<=tiny(asum)) error stop 'bound_calc_green_and_self: zero assignment response'
               green_out(i,j,k)=cmplx(real(bsum/(asum*asum),4),0.0,kind=4)
            enddo
         enddo
      enddo
      !$omp end parallel do
      ft_buf=green_out
      call sfftw_execute_dft_c2r(iplan,ft_buf,rho_buf)
      call bound_normalize_grid(rho_buf,n)
      do iz=-2,2
         do iy=-2,2
            do ix=-2,2
               self_out(ix,iy,iz)=rho_buf(modulo(ix,n)+1,modulo(iy,n)+1,modulo(iz,n)+1)
            enddo
         enddo
      enddo
   endsubroutine bound_calc_green_and_self

   subroutine bound_normalize_grid(rho,n)
      real(4),intent(inout) :: rho(:,:,:)
      integer,intent(in) :: n
      rho=rho/real(n**3,4)
   endsubroutine bound_normalize_grid

   !----------------------------------------------------------------------
   ! bound_check: 单个候选 Halo 的自束缚判定与迭代去束缚主子程序
   !
   ! 参数说明：
   !   ip       [in]  : 候选 Halo 成员单链表的链头索引
   !   np_fof   [in]  : 原始 FoF 搜索给出的初始候选成员数
   !   x_mean   [in]  : 保留接口参数；当前实现不使用，中心由保留成员重新计算
   !   v_mean   [in]  : 保留接口参数；当前实现不使用，均速由保留成员重新计算
   !   np_bound [out] : 当前保留成员数；只有 status=0 才表示已确认自洽收敛
   !   x_bound  [out] : status=0/5 时的保留成员质心；status=1 时不计算，调用方不得使用
   !   v_bound  [out] : 最终保留成员的平均速度（标准格点超共形速度）
   !
   ! 入口依赖已准备好的任务与独占槽，不重查 FoF 链；数值非有限或源域越界仍会终止。
   ! bound_check_batch 始终回收 status；单独调用本入口时，调用方也应读取 status。
   !----------------------------------------------------------------------
   ! [整改说明 G4]: 增加可选第 8 参数 status (0=OK, 1=TOO_FEW, 5=MAX_ITER)
    subroutine bound_check(ip, np_fof, x_mean, v_mean, np_bound, x_bound, v_bound, status,slot)
      implicit none
      integer,intent(in) :: ip,np_fof
      real(4),intent(in) :: x_mean(3),v_mean(3)
      integer,intent(out) :: np_bound
      real(4),intent(out) :: x_bound(3),v_bound(3)
       integer,optional,intent(out) :: status
       integer,optional,intent(in) :: slot

      integer s,gid,st_loc,l
      real(4) dummy_rho(1,1,1)
      complex(4) dummy_ft(1,1,1)

      ! 该入口是批执行层调用的单 Halo 求解入口，必须已将本任务元数据装入线程独占槽。
       s=omp_get_thread_num()+1
       if (present(slot)) s=slot
      gid=bound_mem(s)%method_gid
      l=worker_slot(s)
      ! 直接路径不使用 PM 数组；占位实参仅维持内部统一接口，不会读写这些占位数组。
      if (gid == 0) then
         call bound_check_grid(ip,np_fof,x_mean,v_mean,np_bound,x_bound,v_bound, &
            s,gid,dummy_rho,dummy_ft,st_loc)
      else
         call bound_check_grid(ip,np_fof,x_mean,v_mean,np_bound,x_bound,v_bound, &
            s,gid,bound_storage(gid)%rho(:,:,:,l),bound_storage(gid)%ft(:,:,:,l),st_loc)
      endif

      if (present(status)) status = st_loc
   endsubroutine bound_check

   subroutine bound_check_grid(ip, np_fof, x_mean, v_mean, np_bound, x_bound, v_bound, &
      s, gid, rho, ft, status)
      implicit none
      integer,intent(in) :: ip,np_fof,s,gid
      real(4),intent(in) :: x_mean(3),v_mean(3)
      integer,intent(out) :: np_bound
      real(4),intent(out) :: x_bound(3),v_bound(3)
      real(4),intent(inout) :: rho(:,:,:)
      complex(4),intent(inout) :: ft(:,:,:)
      integer,intent(out) :: status

      integer max_passes
      real(8),parameter :: batch_fraction=0.02d0,refresh_fraction=0.20d0

      integer n,i,j,k,ix,iy,iz,pass,deleted,nu,nb,ios,nmin
      integer left,last,root,child,key,chosen_lev,gsel,lsel,st_sel,gid_eff
      real(8) sp_sel
      ! 粒子数值缓冲由固定工作槽持有，批次开始时预分配。
      ! [阶段 2.2 增量去束缚]: 增量状态管理局部变量
      integer :: n_just_del
      ! [C2串行优化一]: 移除确认轮标志 is_confirm_pass。
      ! 增量势轮在"无删除且无非束缚粒子"时即为自洽收敛(本轮 C1 已对当前成员完整重算,
      ! C2 为增量自洽值), 不再额外执行一次完整 C2 确认轮。
       logical :: is_full_c2
       integer :: eff_threads
       real(4),allocatable :: phi_in(:)
      real(4) scale,spacing,anchor(3),origin(3),lo(3),hi(3),span
      real(4) t(3),dv(3),raw,delta_zero_code,radius_r,c2_pref,delta_phi
      real(8) qc(3),vc(3),h_target,h_split,r2,r_sep,delta_depth,dx_vec(3)
      integer n_source,sx,sy,sz,tx,ty,tz,axis,nearest
      real(4) coord,delta,sw,tw
      real(8) acc_self
      logical greater

      np_bound=0; x_bound=0.0; v_bound=0.0
      status=0
      max_passes=max(1, bound_max_passes)

      scale=real(ratio_cs,4)
      nmin=ceiling(dble(np_halo_min))
      n=0
      if (gid >= 1) n=bn(gid)

      ! 每个新 Halo 重置删除计数和数值状态，不能沿用上一任务在同一槽中的增量势。
       n_just_del = 0
       eff_threads=bound_c2_threads
       if (gid<4 .or. np_fof<bound_c2_min_np) eff_threads=1
       allocate(phi_in(np_fof),stat=ios)
       if (ios/=0) error stop 'bound_check_grid: inner potential allocation failed'
       phi_in=0.0

      associate(member=>bound_mem(s)%member,keep=>bound_mem(s)%keep, &
         order=>bound_mem(s)%order,q=>bound_mem(s)%q(:,1:np_fof), &
         self=>bound_mem(s)%self,phi=>bound_mem(s)%phi,energy=>bound_mem(s)%energy, &
         first_tsc=>bound_mem(s)%first_tsc,wt_tsc=>bound_mem(s)%wt_tsc, &
         phi_c1=>bound_mem(s)%phi_c1,phi_c2=>bound_mem(s)%phi_c2, &
         just_del_list=>bound_mem(s)%just_del_list,phi_direct=>bound_mem(s)%phi_direct)
         phi(1:np_fof)=0.0; energy(1:np_fof)=0.0
         phi_direct(1:np_fof)=0d0

          j=bound_member_offset(ip)
          member(1:np_fof)=bound_member_table(j:j+np_fof-1)


         ! 原地堆排序：按粒子稳定 ID 升序重排成员表，使结果与链表遍历顺序解耦
         ! 排序终止条件：last==1 时写回最后一个成员并触发 exit 退出
         if (np_fof > 1) then
            left=np_fof/2+1; last=np_fof
            do
               if (left>1) then
                  left=left-1; key=member(left)
               else
                  key=member(last); member(last)=member(1); last=last-1
                  if (last==1) then
                     member(1)=key
                     exit
                  endif
               endif
               root=left; child=2*left
               do while (child<=last)
                  if (child<last) then
                     if (pid(member(child))<pid(member(child+1))) child=child+1
                  endif
                  if (pid(key)>=pid(member(child))) exit
                  member(root)=member(child); root=child; child=2*root
               enddo
               member(root)=key
            enddo
         endif

          do i=2,np_fof
             if (pid(member(i))==pid(member(i-1))) error stop 'bound_check_grid: duplicate stable ID'
          enddo
          ! 初始假设全体候选成员都处于束缚状态，迭代中逐步剔除
         keep(1:np_fof)=.true.; np_bound=np_fof
         ! 与分类阶段使用同一链头 anchor；成员按稳定 ID 排序后仍不更换参考点。
         anchor=bound_mem(s)%anchor
         do i=1,np_fof
            q(:,i)=(xv(1:3,member(i))-anchor)*scale
         enddo
         ! 坐标有限性属于 FoF 输入契约。

         ! 执行使用任务准备保存的映射参数。
         span=real(bound_mesh(s)%span,4)

         gid_eff=gid
         spacing=real(bound_mesh(s)%spacing,4)

         ! 将连续坐标映射到局部 PM 网格节点，并预计算 TSC 权重及自身势
         ! 保留最终源域断言：映射坐标必须位于 [0, n/2-1)，选档准备须事先保证覆盖。
         if (gid_eff >= 1) then
            origin=bound_mem(s)%origin
            do i=1,np_fof
               t=(q(:,i)-origin)/spacing
               if (.not.all(ieee_is_finite(t))) then
                  print*,'bound_check: non-finite grid coordinate'
                   error stop 1
               endif
               if (any(t<0.0) .or. any(t>=real(n/2-1,4))) then
                  print*,'bound_check: particle falls outside source domain'
                   error stop 1
               endif
               do axis=1,3
                  coord = t(axis)
                  nearest = nint(coord)
                  delta = coord - real(nearest, 4)
                  first_tsc(axis, i) = nearest
                  wt_tsc(0, axis, i) = 0.5 * (0.5 - delta)**2
                  wt_tsc(1, axis, i) = 0.75 - delta**2
                  wt_tsc(2, axis, i) = 0.5 * (0.5 + delta)**2
               enddo

               acc_self = 0d0
               do sz = 0, 2
                  do sy = 0, 2
                     do sx = 0, 2
                        sw = wt_tsc(sx, 1, i) * wt_tsc(sy, 2, i) * wt_tsc(sz, 3, i)
                        do tz = 0, 2
                           do ty = 0, 2
                              do tx = 0, 2
                                 tw = wt_tsc(tx, 1, i) * wt_tsc(ty, 2, i) * wt_tsc(tz, 3, i)
                                 acc_self = acc_self + dble(sw * tw) * &
                                    dble(bound_gk(gid)%self_phi(tx-sx, ty-sy, tz-sz))
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
               self(i) = real(acc_self, 4)
            enddo

            ! [C2串行优化五]: 粒子位置在本 halo 去束缚全程不变,
            ! cell (CSR) 结构每 halo 只构建一次, 后续所有 pass 复用
             if (gid_eff==4) call bound_two_level_build_halo(s,np_fof,q,spacing,origin,real(app,4))
             if (bound_tl_ctx(s)%n_active==0) call bound_build_cells(s, np_fof, q, spacing)
         endif

         ! 外层刷新当前保留集合的势；内层保持这份势，按新均速/质心分批删除。
          do pass=1,max_passes
            if (gid_eff == 0) then
               n_source = count(keep(1:np_fof))
               delta_zero_code = 0.0
               is_full_c2 = .true.
               c2_pref = real(a_bound / (6d0 * pi), 4)
               ! 首轮全部成员均保留；只访问 i<j，同时给两粒子累加负势，不乘 1/2。
               if (pass == 1) then
                  phi_direct = 0d0
                  do i=1,np_fof-1
                     do j=i+1,np_fof
                        dx_vec = dble(q(:,j) - q(:,i))
                        r_sep = sqrt(sum(dx_vec*dx_vec))
                        delta_depth = (a_bound/(6d0*pi)) * bound_depth_continuous(r_sep, dble(app))
                        phi_direct(i) = phi_direct(i) - delta_depth
                        phi_direct(j) = phi_direct(j) - delta_depth
                     enddo
                  enddo
               else
                  ! 仅在外层刷新时扣除上一外层后删除的源；负势贡献移除后，存活粒子的势变浅。
                  do k=1,n_just_del
                     j=just_del_list(k)
                     do i=1,np_fof
                        if (.not.keep(i)) cycle
                        dx_vec = dble(q(:,j) - q(:,i))
                        r_sep = sqrt(sum(dx_vec*dx_vec))
                        phi_direct(i) = phi_direct(i) + (a_bound/(6d0*pi)) * &
                           bound_depth_continuous(r_sep, dble(app))
                     enddo
                  enddo
               endif
               ! 本次删除列表已消费；phi_direct 留给后续外层，内层只读取当前 phi。
               n_just_del = 0
               do i=1,np_fof
                  if (.not.keep(i)) cycle
                  phi(i) = real(phi_direct(i), 4)
               enddo
            else
               ! C2 首轮完整求和，后续按上一轮删除成员增量扣除。
               ! [C2串行优化一]: 移除确认轮触发条件 is_confirm_pass
               if (pass == 1) then
                  is_full_c2 = .true.
               else
                  is_full_c2 = .false.
               endif

               ! (a) TSC 沉积当前保留成员的质量源项（清零全网格，后半域保持零填充）
               rho=0.0
               do i=1,np_fof
                  if (.not.keep(i)) cycle
                  do sz=0,2
                     iz = first_tsc(3, i) + sz
                     do sy=0,2
                        iy = first_tsc(2, i) + sy
                        do sx=0,2
                           ix = first_tsc(1, i) + sx
                           rho(ix, iy, iz) = rho(ix, iy, iz) + wt_tsc(sx, 1, i) * wt_tsc(sy, 2, i) * wt_tsc(sz, 3, i)
                        enddo
                     enddo
                  enddo
               enddo

               ! (b) 前向 FFT -> 频域逐点乘 C1 势最优核 -> 逆向 FFT 回空间域标量势
               call sfftw_execute_dft_r2c(bound_mesh(s)%plan,rho,ft)
               ft=ft*bound_gk(gid)%green_phi
               call sfftw_execute_dft_c2r(bound_mesh(s)%iplan,ft,rho)
               rho=rho/real(n**3, 4)

               ! (c) 标量势二次样条回插 (TSC)，扣除粒子自身引力 self(i)
               n_source = count(keep(1:np_fof))
               radius_r = real(n, 4) / 2.0
               if (n_source > 1) then
                  delta_zero_code = real(n_source - 1, 4) / radius_r
               else
                  delta_zero_code = 0.0
               endif

               do i=1,np_fof
                  if (.not.keep(i)) cycle
                  raw=0.0
                  do sz=0,2
                     iz = first_tsc(3, i) + sz
                     do sy=0,2
                        iy = first_tsc(2, i) + sy
                        do sx=0,2
                           ix = first_tsc(1, i) + sx
                           raw = raw + rho(ix, iy, iz) * wt_tsc(sx, 1, i) * wt_tsc(sy, 2, i) * wt_tsc(sz, 3, i)
                        enddo
                     enddo
                  enddo
                  phi_c1(i)=real((a_bound/(6d0*pi*dble(spacing))), 4) * ((raw - self(i)) - delta_zero_code)
                  if (is_full_c2) then
                     phi_c2(i)=0.0
                  endif
               enddo

                if (gid_eff==4 .and. bound_tl_ctx(s)%n_active>0) then
                   call bound_two_level_run_pm(s,np_fof,keep,spacing,a_bound,phi_in,eff_threads)
                else
                   phi_in=0.0
                endif
                if (gid_eff==4 .and. bound_tl_ctx(s)%n_active>0) then
                   if (is_full_c2) then
                      call bound_two_level_pp_build(s,np_fof,q,spacing,bound_tl_ctx(s)%particle_r)
                      call bound_two_level_pp_full(s,np_fof,q,keep,spacing,real(app,4),a_bound, &
                         bound_tl_ctx(s)%particle_r,phi_c2,eff_threads)
                   else
                      call bound_two_level_pp_update(s,np_fof,q,keep,just_del_list,n_just_del, &
                         spacing,real(app,4),a_bound,bound_tl_ctx(s)%particle_r,phi_c2,eff_threads)
                   endif
                else
                   if (is_full_c2) then
                      if (eff_threads>1) then
                         call bound_c2_calc_cells_omp(s,np_fof,q,keep,spacing,real(app,4),a_bound,phi_c2,eff_threads)
                      else
                         call bound_c2_calc_cells(s,np_fof,q,keep,spacing,real(app,4),a_bound,phi_c2)
                      endif
                   else
                      if (eff_threads>1) then
                         call bound_c2_update_incremental_cells_omp(s,np_fof,q,keep,just_del_list,n_just_del, &
                            spacing,real(app,4),a_bound,phi_c2,eff_threads)
                      else
                         call bound_c2_update_incremental_cells(s,np_fof,q,keep,just_del_list,n_just_del, &
                            spacing,real(app,4),a_bound,phi_c2)
                      endif
                   endif
                endif
                n_just_del=0
            endif

            ! 删除位置的势不再清零，因为此后所有数值访问均由 keep 过滤；仅当前成员有效。
            do i = 1, np_fof
               if (.not.keep(i)) cycle
                if (gid_eff>=1) then
                   if (gid_eff==4 .and. bound_tl_ctx(s)%n_active>0) then
                      phi(i)=phi_c1(i)+phi_in(i)+phi_c2(i)
                   else
                      phi(i)=phi_c1(i)+phi_c2(i)
                   endif
                endif
               if (np_bound == 1 .or. n_source == 1) phi(i) = 0.0
               if (.not.ieee_is_finite(phi(i))) then
                  print*,'bound_check: non-finite potential value'
                   error stop 1
               endif
            enddo

            deleted=0
            ! (d) 内层循环：基于最新势场分批迭代剔除非束缚成员
            do
               ! 在当前保留集合上重算质心与平均速度；采用双精度累加，再用于相对动能。
               qc=0d0; vc=0d0
               do i=1,np_fof
                  if (.not.keep(i)) cycle
                  qc=qc+dble(q(:,i)); vc=vc+dble(xv(4:6,member(i)))
               enddo
               qc=qc/dble(np_bound); vc=vc/dble(np_bound)

               ! 统计比能量 >= 0 的非束缚成员
               nu=0
               do i=1,np_fof
                  if (.not.keep(i)) cycle
                  dv=xv(4:6,member(i))-real(vc,4)+real(b_bound,4)*(q(:,i)-real(qc,4))
                  energy(i)=sum(dv*dv)/2.0+phi(i)
                  if (.not.ieee_is_finite(energy(i))) then
                     print*,'bound_check: non-finite particle energy'
                      error stop 1
                  endif
                  if (energy(i)>=0.0) then
                     nu=nu+1
                     order(nu)=i
                  endif
               enddo

               ! 判据 1: 当前势下已无非束缚粒子
               if (nu==0) then
                  if (deleted==0) then
                     ! 本外层势已经与当前保留集合对应，且内层未再删除成员，故可确认自洽收敛。
                     ! PM 不另跑全量 C2 确认轮；直接路径也不另跑完整求和确认轮。
                     ! 若本轮已经删除过成员，即使 nu=0 也必须进入下一外层更新旧势。
                     ! 注意：这里判断的是势源集合与 keep 一致，不仅是当前所有能量为负。
                     ! PM 的 C1 已重算、C2 已更新；直接路径的势也已刷新，且本轮没有再删除成员。
                     ! 最终中心只在成功收敛或未确认收敛保留结果时输出。
                     call calc_bound_center
                      status = 0  ! 0 = OK (自洽收敛)
                      call bound_two_level_release_halo(s)
                     return  ! 完全收敛！
                  endif
                  exit  ! 曾删除粒子，旧势含多余质量，退出内层重新求势确认
               endif

               ! 原地堆排序：按比能量确定删除优先级；能量相同时优先删除成员局部索引较小的粒子
               ! 排序终止条件：last==1 时排序完成并触发 exit 退出
               ! nu=1 不可进入现有堆排序，否则 last 会减到0；长度1的待删列表无需排序。
               if (nu>1) then
                  left=nu/2+1; last=nu
                  do
                     if (left>1) then
                        left=left-1; key=order(left)
                     else
                        key=order(last); order(last)=order(1); last=last-1
                        if (last==1) then
                           order(1)=key
                           exit
                        endif
                     endif
                     root=left; child=2*left
                     do while (child<=last)
                        if (child<last) then
                           greater=energy(order(child+1))>energy(order(child))
                           if (energy(order(child+1))==energy(order(child))) &
                              greater=order(child+1)<order(child)
                           if (greater) child=child+1
                        endif
                        greater=energy(order(child))>energy(key)
                        if (energy(order(child))==energy(key)) greater=order(child)<key
                        if (.not.greater) exit
                        order(root)=order(child); root=child; child=2*root
                     enddo
                     order(root)=key
                  enddo
               endif

               ! 单批删除当前剩余成员的 2%（至少 1 个）
               nb=min(nu,max(1,ceiling(batch_fraction*dble(np_bound))))
               ! 这里只登记删除，不立即更新势；保持2%批删除与20%外层刷新规则。
               do i=nu,nu-nb+1,-1
                  keep(order(i))=.false.
                  n_just_del = n_just_del + 1
                  just_del_list(n_just_del) = order(i)
               enddo
               np_bound=np_bound-nb; deleted=deleted+nb

               if (np_bound<nmin) then
                   status = 1  ! 1 = TOO_FEW (剩余成员不足有效门槛，不计算中心)
                   call bound_two_level_release_halo(s)
                  return
               endif

               ! 判据 2: 累计删除数达到当前剩余成员的 20%，退出内层并刷新势。
               if (dble(deleted)>=refresh_fraction*dble(np_bound)) exit
            enddo
         enddo

         ! status=5 是未确认收敛，不是正常通过或物理成员不足；保留成员并输出明确错误报告。
         ! [整改说明 G4]: 达到外层上限未确认收敛，输出当前保留质心并返回 MAX_ITER (5)
         call calc_bound_center
         print*,'bound_check: MAX_ITER, candidate =',ip,' retained members =',np_bound
          status = 5  ! 5 = MAX_ITER
          call bound_two_level_release_halo(s)
         return
      endassociate

   contains

      ! associate 构造内的别名不用于内部子程序，改为显式访问宿主槽的成员缓冲。
      ! 改为显式访问宿主 bound_mem(s) 的成员 keep 与 member，保证作用域合规。
      subroutine calc_bound_center
         real(8) qc_loc(3),vc_loc(3)
         integer ii
         if (np_bound>0) then
            qc_loc=0d0; vc_loc=0d0
            do ii=1,np_fof
               if (.not.bound_mem(s)%keep(ii)) cycle
               qc_loc=qc_loc+dble(bound_mem(s)%q(:,ii)); vc_loc=vc_loc+dble(xv(4:6,bound_mem(s)%member(ii)))
            enddo
            x_bound=anchor+real(qc_loc/(dble(np_bound)*dble(scale)),4)
            v_bound=real(vc_loc/dble(np_bound),4)
         endif
      endsubroutine calc_bound_center

   endsubroutine bound_check_grid

   !----------------------------------------------------------------------
   ! bound_classify: 串行准备任务元数据，并按固定组形成从大到小的执行队列
   !
   ! 参数说明：
   !   head        [in]  : 各候选 Halo 的CRS起点数组 (nh)
   !   np_list     [in]  : 各候选初始成员数数组 (nh)，已满足 FoF 门槛
   !   nthreads    [in]  : 当前 CPU 预算，范围 1..fofcore
   !   ng_standard [in]  : 兼容接口参数；本批入口传单节点 ng，选档使用模块中的模拟参数
   !   schedule    [out] : 预选方法、组别、坐标元数据、组内排序队列与容量
   !   weights     [in, opt] : 工作量打分权重 [A_fft, B_grid, C_part, D_sort]
   !----------------------------------------------------------------------
   subroutine bound_classify(head, np_list, nthreads, ng_standard, schedule, weights)
      implicit none
      integer,intent(in) :: head(:),np_list(:),nthreads,ng_standard
      type(bound_schedule_type),intent(out) :: schedule
      real(8),optional,intent(in) :: weights(4)

      integer nh,ns(4),njob(4),cursor(4),active(4),nt(4)
      integer h,i,j,g,p,nactive,ios,left,last,root,child,key,beg,fin,w,s,best,gsel,lsel,st_sel
      real(8) c(4),nn,volume,score,sp_sel
      real(8) work(4),largest(4),gain,best_gain,before,after
      real(4) anchor(3),lo(3),hi(3),q(3)
      logical greater

      schedule=bound_schedule_type()
       if (nthreads<1 .or. nthreads>fofcore) error stop 'bound_classify: invalid worker count'
      nh=size(head)
      if (size(np_list)/=nh .or. nthreads<1 .or. ng_standard<1) then
          print*,'bound_classify: invalid input array sizes or parameters'
          error stop 1
      endif


      ns=bn(1:4)
      ! 网格类别与 spacing 统一由 bound_select_grid_and_level 按初始 N 和真实跨度 D 选择。

      ! 估算工作量四项权重：正反 FFT、网格扫描、粒子扫描、成员排序
      c=[6d0,12d0,64d0,9d0]
      if (present(weights)) c=weights
      if (.not.all(ieee_is_finite(c)) .or. any(c<0d0) .or. maxval(c)==0d0) then
          print*,'bound_classify: invalid weights'
          error stop 1
      endif

      allocate(schedule%group(nh),schedule%method_gid(nh),schedule%level(nh),schedule%spacing(nh), &
         schedule%flag(nh),schedule%span(nh), &
         schedule%cost(nh),schedule%queue(nh),schedule%slot(nthreads,1), &
         schedule%anchor(3,nh),schedule%origin(3,nh),stat=ios)
      if (ios/=0) then
          print*,'bound_classify: schedule allocation failed'
          error stop 1
      endif
      schedule%group=0; schedule%method_gid=0; schedule%level=0; schedule%spacing=0d0
      schedule%flag=0; schedule%span=0.0
      schedule%cost=0d0; schedule%queue=0; schedule%slot=0
      schedule%anchor=0.0; schedule%origin=0.0
      njob=0; work=0d0; largest=0d0

      ! 1. 遍历每条候选链一次，按真实最大轴向跨度分类；坐标必须已在周期边界展开
      do h=1,nh
          beg=bound_member_offset(h)
          p=bound_member_table(beg)
          anchor=xv(1:3,p); lo=0.0; hi=0.0
          schedule%anchor(:,h)=anchor
          do i=1,np_list(h)
             p=bound_member_table(beg+i-1)
             q=(xv(1:3,p)-anchor)*real(ratio_cs,4)
             lo=min(lo,q); hi=max(hi,q)
          enddo
          schedule%span(h)=maxval(hi-lo)
         if (.not.ieee_is_finite(schedule%span(h))) then
             print*,'bound_classify: non-finite span'
             error stop 1
         endif
         ! 候选真实三维展开跨度严格上限 D_max = 15 Mpc 用于 PM 网格覆盖；超限打印完整配置并以
         ! 非零状态终止，不放粗 spacing 或截断成员；直接求和不受此 PM 覆盖上限限制。
         if (np_list(h) > 1000 .and. schedule%span(h) > real(s_15,4)) then
            call bound_grid_reject('bound_classify',head(h),dble(schedule%span(h)), &
               'candidate span exceeds D_max = 15 Mpc')
         endif
         ! 按正式 N+D 规则选档，取其网格类别作为分组
         ! 唯一方法/网格决策；求解阶段仅使用保存结果，不再调用选档器。
         call bound_select_grid_and_level(int(np_list(h),8), dble(schedule%span(h)), gsel, lsel, sp_sel, st_sel)
         if (st_sel /= BOUND_SEL_OK) then
            call bound_grid_reject('bound_classify',head(h),dble(schedule%span(h)), &
               'no enabled grid class or direct combination can cover candidate')
         endif
         if (gsel == 0) then
            schedule%method_gid(h) = 0
            schedule%level(h) = 0
            schedule%spacing(h) = 0d0
            g = 1  ! 直接求和任务归 48 组 (group 1)
            schedule%group(h) = 1
            nn=dble(np_list(h))
            score=nn*(nn-1d0)/2d0
         else
            schedule%method_gid(h) = gsel
            schedule%level(h) = lsel
            ! 按实际执行的单精度 spacing 再确认覆盖，并保存双精度镜像，避免两阶段单位/精度不一致。
            schedule%spacing(h) = dble(real(sp_sel,4))
            schedule%origin(:,h)=lo-real(sp_sel,4)
            if (dble(bn(gsel)/2-3)*schedule%spacing(h)<dble(schedule%span(h))) &
               call bound_grid_reject('bound_classify',head(h),dble(schedule%span(h)), &
               'actual single precision spacing does not cover span')
            schedule%pm_capacity(gsel)=max(schedule%pm_capacity(gsel),np_list(h))
            g = gsel
            schedule%group(h) = g
            nn=dble(np_list(h)); volume=dble(ns(g))**3
            score=c(1)*volume*log(volume)/log(2d0)+c(2)*volume+ &
               c(3)*nn+c(4)*nn*log(max(2d0,nn))/log(2d0)
         endif
         if (.not.ieee_is_finite(score)) then
             print*,'bound_classify: non-finite cost score'
             error stop 1
         endif
         schedule%cost(h)=score
         njob(g)=njob(g)+1; work(g)=work(g)+score
         largest(g)=max(largest(g),score)
         schedule%capacity(g)=max(schedule%capacity(g),np_list(h))
      enddo
      schedule%work=work

      ! cost 仅决定组内领取先后，不改变选档或固定核心数；它不是已标定的耗时预测。
      ! 2. 四条连续任务队列，各自按 cost 从大到小原地堆排序；并列按输入编号
      cursor(1)=1
      do g=2,4
         cursor(g)=cursor(g-1)+njob(g-1)
      enddo
      schedule%first=cursor; schedule%last=cursor+njob-1
      do h=1,nh
         g=schedule%group(h)
         if (g==0) cycle
         schedule%queue(cursor(g))=h; cursor(g)=cursor(g)+1
      enddo
      do g=1,4
         if (njob(g)<2) cycle
         associate(order=>schedule%queue(schedule%first(g):schedule%last(g)))
            left=njob(g)/2+1; last=njob(g)
            do
               if (left>1) then
                  left=left-1; key=order(left)
               else
                  key=order(last); order(last)=order(1); last=last-1
                  if (last==1) then
                     order(1)=key
                     exit
                  endif
               endif
               root=left; child=2*left
               do while(child<=last)
                  if (child<last) then
                     greater=schedule%cost(order(child+1))<schedule%cost(order(child))
                     if (schedule%cost(order(child+1))==schedule%cost(order(child))) &
                        greater=order(child+1)>order(child)
                     if (greater) child=child+1
                  endif
                  greater=schedule%cost(order(child))<schedule%cost(key)
                  if (schedule%cost(order(child))==schedule%cost(key)) greater=order(child)>key
                  if (.not.greater) exit
                  order(root)=order(child); root=child; child=2*root
               enddo
               order(root)=key
            enddo
         endassociate
      enddo

      ! 工作量仅决定组内顺序；worker 不跨组借用。
      schedule%nworker=bound_workers
       schedule%slot(1:nthreads,1)=worker_group(1:nthreads)
      schedule%nwave=1
   endsubroutine bound_classify

   !----------------------------------------------------------------------
   ! 连续软化解析积分势深度：供直接求和及 C2 残差使用，不包含质量/单位归一化。
   !----------------------------------------------------------------------
   ! 多项式采用 Horner 形式，分段定义、系数与边界 (r >= h, e <= 1) 保持不变。
   ! r <= 0 时直接返回解析极限 208/(70h)，避免除零异常；仅在 r >= h 时计算 1/r。
   pure real(8) function bound_depth_continuous(r, h) result(d)
      implicit none
      real(8),intent(in) :: r, h
      real(8) :: e, primitive

      if (r <= 0d0) then
         d = 208d0 / (70d0 * h)
      else if (r >= h) then
         d = 1d0 / r
      else
         e = 2d0 * r / h
         if (e <= 1d0) then
            ! Horner: e^2*(112 + e^2*(-56 + e*(14 + e*(8 - 3e))))
            primitive = e*e * (112d0 + e*e*(-56d0 + e*(14d0 + e*(8d0 - 3d0*e))))
            d = (208d0 - primitive) / (70d0 * h)
         else
            ! Horner: e^2*(448 + e*(-280 + e*(56 + e*(14 + e*(-8 + e))))) - 224e - 12/e
            primitive = e*e * (448d0 + e*(-280d0 + e*(56d0 + e*(14d0 + e*(-8d0 + e))))) &
               - 224d0*e - 12d0/e
            d = (128d0 - primitive) / (70d0 * h)
         endif
      endif
   endfunction bound_depth_continuous

   ! C2 软化势残差：计算目标软化势与分裂软化势深度之差。
   ! 此处 r 与 h 均以局部 PM 格长为单位；直接求和则使用标准格坐标和 app。
   pure real(8) function bound_c2_delta_depth(r, h_target, h_split) result(dd)
      implicit none
      real(8),intent(in) :: r, h_target, h_split

      if (r <= 0d0) then
         dd = (208d0 / 70d0) * (1d0 / h_target - 1d0 / h_split)
      else if (r >= h_target .and. r >= h_split) then
         dd = 0d0
      else
         dd = bound_depth_continuous(r, h_target) - bound_depth_continuous(r, h_split)
      endif
   endfunction bound_c2_delta_depth

   !----------------------------------------------------------------------
   ! [C2串行优化五]: 构建 每 halo 一次 的连续 cell (CSR) 结构
   ! u = q/spacing, cell 宽 cw = 3.5 (与链表版一致); 覆盖全部 np_fof 粒子
   ! (保证删除粒子的 cell_id 在增量阶段始终有效)。
   ! 计数排序 O(N) 填充，cell 内按成员局部索引升序；不足时仅扩展本槽容量。
   !----------------------------------------------------------------------
   subroutine bound_build_cells(s, np_fof, q, spacing)
      implicit none
      integer,intent(in) :: s, np_fof
      real(4),intent(in) :: q(3, np_fof), spacing
      real(8),parameter :: cw = 3.5d0
      integer :: i, c, nt, ic(3), ncell(3), need
       integer :: fill_ptr,ntask,p_cnt,p_cur,p_next
      real(8) :: u(3), u_min(3), u_max(3)

      u_min = 1d30; u_max = -1d30
      do i = 1, np_fof
         u = dble(q(:, i)) / dble(spacing)
         u_min = min(u_min, u); u_max = max(u_max, u)
      enddo
      do i = 1, 3
         ncell(i) = max(1, int(floor((u_max(i) - u_min(i)) / cw)) + 1)
      enddo
      nt = ncell(1) * ncell(2) * ncell(3)

      ! [C2审核修正] bound_cells 外层为编译期静态数组, 此处不再有任何共享 allocate
      need = max(np_fof, 1)
      if (.not.allocated(bound_cells(s)%cell_id)) then
         allocate(bound_cells(s)%cell_id(0),bound_cells(s)%cell_particles(0))
      endif
      if (size(bound_cells(s)%cell_id) < need) then
         if (allocated(bound_cells(s)%cell_id)) deallocate(bound_cells(s)%cell_id)
         if (allocated(bound_cells(s)%cell_particles)) deallocate(bound_cells(s)%cell_particles)
         allocate(bound_cells(s)%cell_id(need), bound_cells(s)%cell_particles(need))
      endif
      if (.not.allocated(bound_cells(s)%cell_start)) allocate(bound_cells(s)%cell_start(0))
      if (size(bound_cells(s)%cell_start) < nt + 1) then
         if (allocated(bound_cells(s)%cell_start)) deallocate(bound_cells(s)%cell_start)
         allocate(bound_cells(s)%cell_start(nt + 1))
      endif

      bound_cells(s)%ncell = ncell
      bound_cells(s)%ncell_total = nt
      bound_cells(s)%u_min = u_min
      bound_cells(s)%cell_start(1:nt + 1) = 0

      ! 计数: 每粒子归属 cell (位置不变, 与 keep 状态无关)
      do i = 1, np_fof
         u = dble(q(:, i)) / dble(spacing)
         ic = min(max(int(floor((u - u_min) / cw)) + 1, 1), ncell)
         c = ic(1) + (ic(2) - 1) * ncell(1) + (ic(3) - 1) * ncell(1) * ncell(2)
         bound_cells(s)%cell_id(i) = c
         bound_cells(s)%cell_start(c + 1) = bound_cells(s)%cell_start(c + 1) + 1
      enddo
      ! 前缀和: cell_start(c+1) 暂为 exclusive 结束位置
      do c = 1, nt
         bound_cells(s)%cell_start(c + 1) = bound_cells(s)%cell_start(c + 1) + bound_cells(s)%cell_start(c)
      enddo
      ! [修正] prefix -> 1-based 起始游标: cell_start(1)=1, cell_start(c)=prefix(c-1)+1
      ! (缺少此步会使 cell 1 从位置 0 越界写入且全体区间错位 1)
      bound_cells(s)%cell_start(1) = 1
      do c = 2, nt + 1
         bound_cells(s)%cell_start(c) = bound_cells(s)%cell_start(c) + 1
      enddo
      ! 填充 (i 升序 -> cell 内升序): cell_start(c) 充当写指针
      do i = 1, np_fof
         c = bound_cells(s)%cell_id(i)
         fill_ptr = bound_cells(s)%cell_start(c)
         bound_cells(s)%cell_particles(fill_ptr) = i
         bound_cells(s)%cell_start(c) = fill_ptr + 1
      enddo
      ! 填充后 cell_start(c) = start(c)+cnt(c) = start(c+1), 回移一位重建 CSR 起始偏移
      do c = nt, 2, -1
         bound_cells(s)%cell_start(c) = bound_cells(s)%cell_start(c - 1)
      enddo
       bound_cells(s)%cell_start(1) = 1
       ntask=0
       do c=1,nt
          p_cnt=bound_cells(s)%cell_start(c+1)-bound_cells(s)%cell_start(c)
          if (p_cnt<=0) cycle
          ntask=ntask+(p_cnt+1023)/1024
       enddo
       if (.not.allocated(bound_cells(s)%task_c1)) then
          allocate(bound_cells(s)%task_c1(0),bound_cells(s)%task_p1(0),bound_cells(s)%task_p2(0))
       endif
       if (size(bound_cells(s)%task_c1)<ntask) then
          deallocate(bound_cells(s)%task_c1,bound_cells(s)%task_p1,bound_cells(s)%task_p2)
          allocate(bound_cells(s)%task_c1(ntask),bound_cells(s)%task_p1(ntask),bound_cells(s)%task_p2(ntask))
       endif
       bound_cells(s)%ntasks=ntask
       ntask=0
       do c=1,nt
          p_cnt=bound_cells(s)%cell_start(c+1)-bound_cells(s)%cell_start(c)
          if (p_cnt<=0) cycle
          p_cur=bound_cells(s)%cell_start(c)
          do while (p_cur<=bound_cells(s)%cell_start(c+1)-1)
             p_next=min(p_cur+1023,bound_cells(s)%cell_start(c+1)-1)
             ntask=ntask+1
             bound_cells(s)%task_c1(ntask)=c
             bound_cells(s)%task_p1(ntask)=p_cur
             bound_cells(s)%task_p2(ntask)=p_next
             p_cur=p_next+1
          enddo
       enddo
    endsubroutine bound_build_cells

   subroutine bound_c2_calc_cells_omp(s,np_fof,q,keep,spacing,app_val,a_val,phi_c2,nth)
      integer,intent(in) :: s,np_fof,nth
      real(4),intent(in) :: q(3,np_fof),spacing,app_val
      logical,intent(in) :: keep(np_fof)
      real(8),intent(in) :: a_val
      real(4),intent(inout) :: phi_c2(np_fof)
      real(8),parameter :: r2_cut=12.25d0
      integer :: n_k,i,p1,p2,i1,i2,c1,c2i,sz,cx,cy,cz,nx,ny,nz,nc(3),tid,tid_k,itask,p1_beg,p1_end
      real(8) :: dx,dy,dz,h_target,h_split,r2,r_sep,delta_depth,r2_cut_raw
      real(4) :: c2_pref,delta_phi,phi_sum
      phi_c2=0.0
      n_k=count(keep(1:np_fof))
      if (n_k<=1 .or. real(spacing)<=0.0) return
      h_target=dble(app_val)/dble(spacing); h_split=bound_h_split
      r2_cut_raw=r2_cut*dble(spacing)*dble(spacing)
      c2_pref=real(a_val/(6d0*pi*dble(spacing)),4)
      nc=bound_cells(s)%ncell
      bound_mem(s)%phi_c2_team(1:np_fof,1:nth)=0.0
      !$omp parallel num_threads(nth) default(shared) &
      !$omp private(tid,tid_k,itask,p1_beg,p1_end,c1,cx,cy,cz,p1,i1,p2,i2,dx,dy,dz, &
      !$omp r2,r_sep,delta_depth,delta_phi,sz,nx,ny,nz,c2i,i,phi_sum)
      !$omp single
      if (omp_get_num_threads()/=nth) error stop 'bound_c2_calc_cells_omp: actual team size mismatch'
      !$omp end single
      tid=omp_get_thread_num()+1
      !$omp do schedule(dynamic,1)
      do itask=1,bound_cells(s)%ntasks
         c1=bound_cells(s)%task_c1(itask)
         p1_beg=bound_cells(s)%task_p1(itask)
         p1_end=bound_cells(s)%task_p2(itask)
         cz=(c1-1)/(nc(1)*nc(2))+1
         cy=mod(c1-1,nc(1)*nc(2))/nc(1)+1
         cx=mod(mod(c1-1,nc(1)*nc(2)),nc(1))+1
         do p1=p1_beg,p1_end
            i1=bound_cells(s)%cell_particles(p1)
            if (.not.keep(i1)) cycle
            do p2=p1+1,bound_cells(s)%cell_start(c1+1)-1
               i2=bound_cells(s)%cell_particles(p2)
               if (.not.keep(i2)) cycle
               dx=dble(q(1,i2)-q(1,i1))
               dy=dble(q(2,i2)-q(2,i1))
               dz=dble(q(3,i2)-q(3,i1))
               r2=dx*dx+dy*dy+dz*dz
               if (r2<r2_cut_raw) then
                  r_sep=sqrt(r2)/dble(spacing)
                  delta_depth=bound_c2_delta_depth(r_sep,h_target,h_split)
                  delta_phi=-c2_pref*real(delta_depth,4)
                  bound_mem(s)%phi_c2_team(i1,tid)=bound_mem(s)%phi_c2_team(i1,tid)+delta_phi
                  bound_mem(s)%phi_c2_team(i2,tid)=bound_mem(s)%phi_c2_team(i2,tid)+delta_phi
               endif
            enddo
            do sz=1,13
               nx=cx+bound_c2_stencil(1,sz)
               ny=cy+bound_c2_stencil(2,sz)
               nz=cz+bound_c2_stencil(3,sz)
               if (nx<1 .or. nx>nc(1) .or. ny<1 .or. ny>nc(2) .or. nz<1 .or. nz>nc(3)) cycle
               c2i=nx+(ny-1)*nc(1)+(nz-1)*nc(1)*nc(2)
               if (bound_cells(s)%cell_start(c2i)>bound_cells(s)%cell_start(c2i+1)-1) cycle
               do p2=bound_cells(s)%cell_start(c2i),bound_cells(s)%cell_start(c2i+1)-1
                  i2=bound_cells(s)%cell_particles(p2)
                  if (.not.keep(i2)) cycle
                  dx=dble(q(1,i2)-q(1,i1))
                  dy=dble(q(2,i2)-q(2,i1))
                  dz=dble(q(3,i2)-q(3,i1))
                  r2=dx*dx+dy*dy+dz*dz
                  if (r2<r2_cut_raw) then
                     r_sep=sqrt(r2)/dble(spacing)
                     delta_depth=bound_c2_delta_depth(r_sep,h_target,h_split)
                     delta_phi=-c2_pref*real(delta_depth,4)
                     bound_mem(s)%phi_c2_team(i1,tid)=bound_mem(s)%phi_c2_team(i1,tid)+delta_phi
                     bound_mem(s)%phi_c2_team(i2,tid)=bound_mem(s)%phi_c2_team(i2,tid)+delta_phi
                  endif
               enddo
            enddo
         enddo
      enddo
      !$omp end do
      !$omp do schedule(static)
      do i=1,np_fof
         phi_sum=0.0
         do tid_k=1,nth
            phi_sum=phi_sum+bound_mem(s)%phi_c2_team(i,tid_k)
         enddo
         phi_c2(i)=phi_sum
      enddo
      !$omp end do
      !$omp end parallel
   endsubroutine bound_c2_calc_cells_omp

   subroutine bound_c2_update_incremental_cells_omp(s,np_fof,q,keep,just_del_list,n_just_del, &
      spacing,app_val,a_val,phi_c2,nth)
      integer,intent(in) :: s,np_fof,n_just_del,nth,just_del_list(n_just_del)
      real(4),intent(in) :: q(3,np_fof),spacing,app_val
      logical,intent(in) :: keep(np_fof)
      real(8),intent(in) :: a_val
      real(4),intent(inout) :: phi_c2(np_fof)
      real(8),parameter :: r2_cut=12.25d0
      integer :: n_k,k,p,j,i2,cx,cy,cz,nx,ny,nz,c_j,c2i,nc(3),tid,tid_k,i
      real(8) :: ddx,ddy,ddz,h_target,h_split,r2,r_sep,delta_depth,r2_cut_raw
      real(4) :: c2_pref,delta_phi,phi_sum
      n_k=count(keep(1:np_fof))
      if (n_just_del<=0 .or. n_k<=0 .or. real(spacing)<=0.0) return
      h_target=dble(app_val)/dble(spacing); h_split=bound_h_split
      r2_cut_raw=r2_cut*dble(spacing)*dble(spacing)
      c2_pref=real(a_val/(6d0*pi*dble(spacing)),4)
      nc=bound_cells(s)%ncell
      bound_mem(s)%phi_c2_team(1:np_fof,1:nth)=0.0
      !$omp parallel num_threads(nth) default(shared) &
      !$omp private(tid,tid_k,k,j,c_j,cz,cy,cx,nz,ny,nx,c2i,p,i2,ddx,ddy,ddz,r2,r_sep,delta_depth,delta_phi,i,phi_sum)
      !$omp single
      if (omp_get_num_threads()/=nth) error stop 'bound_c2_update_cells_omp: actual team size mismatch'
      !$omp end single
      tid=omp_get_thread_num()+1
      !$omp do schedule(dynamic,4)
      do k=1,n_just_del
         j=just_del_list(k)
         c_j=bound_cells(s)%cell_id(j)
         cz=(c_j-1)/(nc(1)*nc(2))+1
         cy=mod(c_j-1,nc(1)*nc(2))/nc(1)+1
         cx=mod(mod(c_j-1,nc(1)*nc(2)),nc(1))+1
         do nz=cz-1,cz+1
            if (nz<1 .or. nz>nc(3)) cycle
            do ny=cy-1,cy+1
               if (ny<1 .or. ny>nc(2)) cycle
               do nx=cx-1,cx+1
                  if (nx<1 .or. nx>nc(1)) cycle
                  c2i=nx+(ny-1)*nc(1)+(nz-1)*nc(1)*nc(2)
                  do p=bound_cells(s)%cell_start(c2i),bound_cells(s)%cell_start(c2i+1)-1
                     i2=bound_cells(s)%cell_particles(p)
                     if (.not.keep(i2)) cycle
                     ddx=dble(q(1,i2)-q(1,j))
                     ddy=dble(q(2,i2)-q(2,j))
                     ddz=dble(q(3,i2)-q(3,j))
                     r2=ddx*ddx+ddy*ddy+ddz*ddz
                     if (r2<r2_cut_raw) then
                        r_sep=sqrt(r2)/dble(spacing)
                        delta_depth=bound_c2_delta_depth(r_sep,h_target,h_split)
                        delta_phi=-c2_pref*real(delta_depth,4)
                        bound_mem(s)%phi_c2_team(i2,tid)=bound_mem(s)%phi_c2_team(i2,tid)-delta_phi
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end do
      !$omp do schedule(static)
      do i=1,np_fof
         if (.not.keep(i)) cycle
         phi_sum=0.0
         do tid_k=1,nth
            phi_sum=phi_sum+bound_mem(s)%phi_c2_team(i,tid_k)
         enddo
         if (phi_sum/=0.0) phi_c2(i)=phi_c2(i)+phi_sum
      enddo
      !$omp end do
      !$omp end parallel
   endsubroutine bound_c2_update_incremental_cells_omp
   ! [C2串行优化二]: 截断平方改在原始 q 坐标域比较 (r2_cut_raw = 12.25 * spacing^2,
   ! 每轮预计算一次)，候选对循环不再做逐对除法；
   ! 仅严格进入截断范围的真实近邻对才计算一次无量纲距离 r = sqrt(r2_raw) / spacing。
   ! 严格 r < 3.5 判据不变：数学上 r_raw < 3.5*spacing 与 (r/spacing)^2 < 12.25 等价。
   ! 浮点舍入路径变化可能影响截断边界附近的对；覆盖性与数值差异需由测试核验。

   ! [C2串行优化三]: 显式标量化, 避免长度3数组临时量, 改善寄存器使用;
   ! 精度顺序保持: 单精度坐标差后转双精度, 平方和按 dx^2+dy^2+dz^2 顺序累加

   !----------------------------------------------------------------------
   ! [C2串行优化五]: 连续 cell (CSR) 版全量 C2 短程解析势
   ! 遍历本 cell 内部 + 13 个半邻域 cell，跳过 keep=false 粒子；
   ! cell 内按成员局部索引升序，每对只处理一次，不保存邻居对列表。
   ! cell 宽 cw=3.5 与截断半径相同，因此所有 r<3.5 的对只能位于本 cell 或相邻 cell。
   ! 13 个半邻域与本 cell 内遍历共同覆盖每个 r<3.5 的粒子对，且只访问一次。
   !----------------------------------------------------------------------
   subroutine bound_c2_calc_cells(s, np_fof, q, keep, spacing, app_val, a_val, phi_c2)
      implicit none
      integer,intent(in) :: s, np_fof
      real(4),intent(in) :: q(3, np_fof)
      logical,intent(in) :: keep(np_fof)
      real(4),intent(in) :: spacing, app_val
      real(8),intent(in) :: a_val
      real(4),intent(inout) :: phi_c2(np_fof)

      real(8),parameter :: r2_cut = 12.25d0
      integer :: n_k, i, p1, p2, i1, i2, c1, c2i, sz, cx, cy, cz, nx, ny, nz
      integer :: nc(3)
      real(8) :: dx, dy, dz, h_target, h_split, r2, r_sep, delta_depth, r2_cut_raw
      real(4) :: c2_pref, delta_phi

      phi_c2 = 0.0
      n_k = count(keep(1:np_fof))

      if (n_k <= 1 .or. real(spacing) <= 0.0) then
         return
      endif

      h_target = dble(app_val) / dble(spacing); h_split = bound_h_split
      r2_cut_raw = r2_cut * dble(spacing) * dble(spacing)
      c2_pref = real(a_val / (6d0 * pi * dble(spacing)), 4)
      nc = bound_cells(s)%ncell(1:3)

      do cz = 1, nc(3)
         do cy = 1, nc(2)
            do cx = 1, nc(1)
               c1 = cx + (cy - 1) * nc(1) + (cz - 1) * nc(1) * nc(2)
               if (bound_cells(s)%cell_start(c1) > bound_cells(s)%cell_start(c1 + 1) - 1) cycle
               do p1 = bound_cells(s)%cell_start(c1), bound_cells(s)%cell_start(c1 + 1) - 1
                  i1 = bound_cells(s)%cell_particles(p1)
                  if (.not. keep(i1)) cycle
                  ! 本 cell 内部: p2 从 p1 的下一位置起 (与链表版一致, 不含 self 对), 每对一次
                  do p2 = p1 + 1, bound_cells(s)%cell_start(c1 + 1) - 1
                     i2 = bound_cells(s)%cell_particles(p2)
                     if (.not. keep(i2)) cycle
                     dx = dble(q(1, i2) - q(1, i1))
                     dy = dble(q(2, i2) - q(2, i1))
                     dz = dble(q(3, i2) - q(3, i1))
                     r2 = dx*dx + dy*dy + dz*dz
                     if (r2 < r2_cut_raw) then
                        r_sep = sqrt(r2) / dble(spacing)
                        delta_depth = bound_c2_delta_depth(r_sep, h_target, h_split)
                        delta_phi = -c2_pref * real(delta_depth, 4)
                        phi_c2(i1) = phi_c2(i1) + delta_phi
                        phi_c2(i2) = phi_c2(i2) + delta_phi
                     endif
                  enddo
                  ! 13 个半邻域 cell
                  do sz = 1, 13
                     nx = cx + bound_c2_stencil(1, sz)
                     ny = cy + bound_c2_stencil(2, sz)
                     nz = cz + bound_c2_stencil(3, sz)
                     if (nx < 1 .or. nx > nc(1) .or. ny < 1 .or. ny > nc(2) .or. nz < 1 .or. nz > nc(3)) cycle
                     c2i = nx + (ny - 1) * nc(1) + (nz - 1) * nc(1) * nc(2)
                     if (bound_cells(s)%cell_start(c2i) > bound_cells(s)%cell_start(c2i + 1) - 1) cycle
                     do p2 = bound_cells(s)%cell_start(c2i), bound_cells(s)%cell_start(c2i + 1) - 1
                        i2 = bound_cells(s)%cell_particles(p2)
                        if (.not. keep(i2)) cycle
                        dx = dble(q(1, i2) - q(1, i1))
                        dy = dble(q(2, i2) - q(2, i1))
                        dz = dble(q(3, i2) - q(3, i1))
                        r2 = dx*dx + dy*dy + dz*dz
                        if (r2 < r2_cut_raw) then
                           r_sep = sqrt(r2) / dble(spacing)
                           delta_depth = bound_c2_delta_depth(r_sep, h_target, h_split)
                           delta_phi = -c2_pref * real(delta_depth, 4)
                           phi_c2(i1) = phi_c2(i1) + delta_phi
                           phi_c2(i2) = phi_c2(i2) + delta_phi
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

   endsubroutine bound_c2_calc_cells

   !----------------------------------------------------------------------
   ! [C2串行优化五]: 连续 cell (CSR) 版增量 C2 更新
   ! 按删除粒子的 cell_id 遍历周围 27 邻格，只更新仍保留的粒子；
   ! 每个 del-keep 交叉对恰好处理一次，不计算 del-del 与 keep-keep。
   !----------------------------------------------------------------------
   subroutine bound_c2_update_incremental_cells(s, np_fof, q, keep, just_del_list, n_just_del, &
      spacing, app_val, a_val, phi_c2)
      implicit none
      integer,intent(in) :: s, np_fof, n_just_del
      real(4),intent(in) :: q(3, np_fof)
      logical,intent(in) :: keep(np_fof)
      integer,intent(in) :: just_del_list(n_just_del)
      real(4),intent(in) :: spacing, app_val
      real(8),intent(in) :: a_val
      real(4),intent(inout) :: phi_c2(np_fof)

      real(8),parameter :: r2_cut = 12.25d0
      integer :: n_k, k, p, j, i2, cx, cy, cz, nx, ny, nz, c_j, c2i
      integer :: nc(3)
      real(8) :: ddx, ddy, ddz, h_target, h_split, r2, r_sep, delta_depth, r2_cut_raw
      real(4) :: c2_pref, delta_phi
      n_k = count(keep(1:np_fof))

      if (n_just_del <= 0 .or. n_k <= 0 .or. real(spacing) <= 0.0) then
         return
      endif

      h_target = dble(app_val) / dble(spacing); h_split = bound_h_split
      r2_cut_raw = r2_cut * dble(spacing) * dble(spacing)
      c2_pref = real(a_val / (6d0 * pi * dble(spacing)), 4)
      nc = bound_cells(s)%ncell(1:3)

      do k = 1, n_just_del
         j = just_del_list(k)
         c_j = bound_cells(s)%cell_id(j)
         ! 1D cell 编号还原 3D 索引
         cz = (c_j - 1) / (nc(1) * nc(2)) + 1
         cy = mod(c_j - 1, nc(1) * nc(2)) / nc(1) + 1
         cx = mod(mod(c_j - 1, nc(1) * nc(2)), nc(1)) + 1
         ! 遍历中心 cell 及周围 26 邻格 (共 27)
         do nz = cz - 1, cz + 1
            if (nz < 1 .or. nz > nc(3)) cycle
            do ny = cy - 1, cy + 1
               if (ny < 1 .or. ny > nc(2)) cycle
               do nx = cx - 1, cx + 1
                  if (nx < 1 .or. nx > nc(1)) cycle
                  c2i = nx + (ny - 1) * nc(1) + (nz - 1) * nc(1) * nc(2)
                  do p = bound_cells(s)%cell_start(c2i), bound_cells(s)%cell_start(c2i + 1) - 1
                     i2 = bound_cells(s)%cell_particles(p)
                     if (.not. keep(i2)) cycle
                     ddx = dble(q(1, i2) - q(1, j))
                     ddy = dble(q(2, i2) - q(2, j))
                     ddz = dble(q(3, i2) - q(3, j))
                     r2 = ddx*ddx + ddy*ddy + ddz*ddz
                     if (r2 < r2_cut_raw) then
                        r_sep = sqrt(r2) / dble(spacing)
                        delta_depth = bound_c2_delta_depth(r_sep, h_target, h_split)
                        delta_phi = -c2_pref * real(delta_depth, 4)
                        ! 扣除删除粒子 j 对存活粒子 i2 的引力吸引贡献 (势变浅)
                        phi_c2(i2) = phi_c2(i2) - delta_phi
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo

   endsubroutine bound_c2_update_incremental_cells
   ! [C2串行优化二]: 截断平方在原始 q 坐标域比较 (r2_cut_raw = 12.25 * spacing^2,
   ! 每轮预计算一次)，候选对循环不再做逐对除法；
   ! 仅严格进入截断范围的真实近邻对才计算无量纲距离 r = sqrt(r2_raw) / spacing。
   ! 严格 r < 3.5 判据不变：数学上 r_raw < 3.5*spacing 与 (r/spacing)^2 < 12.25 等价。
   ! 浮点舍入路径变化可能影响截断边界附近的对；覆盖性与数值差异需由测试核验。

   ! 确定 cell 范围 (覆盖全部 np_fof 粒子，保证任意删除粒子与剩余粒子的 cell 坐标有效落在 [1, ncell])

! [C2串行优化三]: 显式标量化 (同全量路径)

   pure real(8) function bound_fft_wave_number(index_1based, ngrid) result(kval)
      implicit none
      integer,intent(in) :: index_1based, ngrid
      integer :: mode
      mode = modulo(index_1based - 1 + ngrid/2, ngrid) - ngrid/2
      kval = 2d0*pi*dble(mode)/dble(ngrid)
   endfunction bound_fft_wave_number

   pure real(8) function bound_assignment_window_sq(kvec, assignment_order) result(value)
      implicit none
      real(8),intent(in) :: kvec(3)
      integer,intent(in) :: assignment_order
      integer :: axis
      real(8) :: x, sinc_product
      sinc_product = 1d0
      do axis = 1, 3
         x = 0.5d0*kvec(axis)
         if (x /= 0d0) sinc_product = sinc_product*sin(x)/x
      enddo
      value = sinc_product**(2*assignment_order + 2)
   endfunction bound_assignment_window_sq

endmodule
