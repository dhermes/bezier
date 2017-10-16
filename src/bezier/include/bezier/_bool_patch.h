// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef BEZIER_BOOL_PATCH_H
#define BEZIER_BOOL_PATCH_H

#if defined (__cplusplus)
extern "C" {
#elif !defined (_MSC_VER)
#include <stdbool.h>
#endif

#if !defined (__cplusplus) && defined (_MSC_VER)
#  if !defined (FALSE)
#    define FALSE (0)
#  endif
#  if !defined (TRUE)
#    define TRUE (!FALSE)
#  endif
#  if _MSC_VER >= 1800
#    include <stdbool.h>
#  else
//   NOTE: We must use `char` because it is 1 byte, as `bool` is
//         in C99 and later (which is what `gfortran` uses).
#    define bool char
#    define true TRUE
#    define false FALSE
#  endif
#endif

#if defined (__cplusplus)
}
#endif

#endif /* BEZIER_BOOL_PATCH_H */
