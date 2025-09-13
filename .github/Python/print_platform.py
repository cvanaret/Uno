import platform
import sysconfig

print(platform.machine()) # x86_64 / i386 / arm64 / AMD64
print(platform.system()) # linux / Darwin / Windows
print(sysconfig.get_platform()) # linux-x86_64 / macosx-10.13-universal2 / win-amd64
print(sys.version) # 3.12.3 (main, Aug 14 2025, 17:47:21) [GCC 13.3.0]