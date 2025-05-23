import os
import sys
import subprocess
import pathlib
import shutil
from setuptools import setup, find_packages
from setuptools.command.build_py import build_py

ext = "dylib" if sys.platform == "darwin" else "so"
libopentrustregion_file = f"libopentrustregion.{ext}"
libtestsuite_file = f"libtestsuite.{ext}"
package_dir = pathlib.Path(__file__).parent.absolute()
build_dir = package_dir / "build"


class CMakeBuild(build_py):
    def run(self):
        # ensure build directory exists
        os.makedirs(build_dir, exist_ok=True)

        # run CMake configure & build
        subprocess.check_call(["cmake", ".."], cwd=build_dir)
        subprocess.check_call(["cmake", "--build", "."], cwd=build_dir)

        # copy built binaries to package directory
        shutil.copy(
            build_dir / libopentrustregion_file,
            package_dir / "pyopentrustregion" / libopentrustregion_file,
        )
        shutil.copy(
            build_dir / libtestsuite_file,
            package_dir / "pyopentrustregion" / libtestsuite_file,
        )

        # run steps in parent class
        super().run()


setup(
    packages=find_packages(),
    include_package_data=True,
    package_data={"pyopentrustregion": [libopentrustregion_file, libtestsuite_file]},
    cmdclass={"build_py": CMakeBuild},
)
