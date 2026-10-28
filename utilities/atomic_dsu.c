#include <stdint.h>
/* 完整复制自 fof_new/clean/atomic_dsu.c；原子操作与根选择不变。 */

/*
 * Lock-free DSU for concurrent FoF halo grouping.
 * Roots are marked by negative values (-1).
 * Non-root entries point to their 1-based parent index.
 * Smaller index is chosen as root ("大挂小").
 */
static int32_t root(int32_t *parent, int32_t p)
{
   for (;;) {
      int32_t next = __atomic_load_n(&parent[p - 1], __ATOMIC_RELAXED);
      if (next < 0) return p;
      int32_t after = __atomic_load_n(&parent[next - 1], __ATOMIC_RELAXED);
      if (after > 0) {
         int32_t expected = next;
         (void)__atomic_compare_exchange_n(&parent[p - 1], &expected, after, 0,
                                           __ATOMIC_RELAXED, __ATOMIC_RELAXED);
      }
      p = next;
   }
}

void fof_atomic_union(int32_t *parent, int32_t i, int32_t j)
{
   for (;;) {
      int32_t a = root(parent, i);
      int32_t b = root(parent, j);
      if (a == b) return;
      int32_t low = a < b ? a : b;
      int32_t high = a < b ? b : a;
      int32_t expected = -1;
      if (__atomic_compare_exchange_n(&parent[high - 1], &expected, low, 0,
                                      __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {
         return;
      }
   }
}
