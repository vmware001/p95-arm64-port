/* Minimal stubs for MSVC runtime symbols referenced by the official prebuilt
 * gwnum64.lib, so the archive can be linked by MinGW for benchmarking.
 * Only used by the benchmark harness (never shipped). */
#include <stdint.h>

extern "C" {
uintptr_t __security_cookie = 0x00002B992DDFA232ull;

void __security_check_cookie(uintptr_t cookie) { (void) cookie; }
void __GSHandlerCheck(void) { }
}
