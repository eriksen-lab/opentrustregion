import os
import sys
import subprocess
import pathlib
import shutil
from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Optional

ext = "dylib" if sys.platform == "darwin" else "so"
package_dir = pathlib.Path(__file__).parent.absolute()
build_dir = package_dir / "build"
libopentrustregion_file: Optional[str]
for suffix in ["", "32", "64"]:
    libopentrustregion_file = f"libopentrustregion{suffix}.{ext}"
    libopentrustregion_path = build_dir / libopentrustregion_file
    if libopentrustregion_path.exists():
        break
else:
    libopentrustregion_file = None
libtestsuite_file = f"libtestsuite.{ext}"


class CMakeBuild(build_py):
    def run(self):
        # ensure build directory exists
        os.makedirs(build_dir, exist_ok=True)

        # add extra flags
        extra_flags = os.getenv("CMAKE_FLAGS", "")
        cmake_cmd = ["cmake", ".."]
        if extra_flags:
            cmake_cmd.append(extra_flags)

        # run CMake configure & build
        subprocess.check_call(cmake_cmd, cwd=build_dir)
        subprocess.check_call(["cmake", "--build", "."], cwd=build_dir)

        # copy libopentrustregion only if it exists (i.e., shared build)
        if libopentrustregion_file is not None:
            shutil.copy(
                libopentrustregion_path,
                package_dir / "pyopentrustregion" / libopentrustregion_file,
            )

        # always copy testsuite
        shutil.copy(
            build_dir / libtestsuite_file,
            package_dir / "pyopentrustregion" / libtestsuite_file,
        )

        # run steps in parent class
        super().run()


# update package_data dynamically
package_data_files = [libtestsuite_file]
if libopentrustregion_file is not None:
    package_data_files.append(libopentrustregion_file)

setup(
    packages=find_packages(),
    include_package_data=True,
    package_data={"pyopentrustregion": package_data_files},
    cmdclass={"build_py": CMakeBuild},
)
