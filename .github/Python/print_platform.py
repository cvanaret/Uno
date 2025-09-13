import platform
import sysconfig

print(platform.machine()) # x86_64 / i386 / arm64 / AMD64
print(platform.system()) # linux / Darwin
print(sysconfig.get_platform()) # linux-x86_64