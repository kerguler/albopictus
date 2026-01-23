from __future__ import annotations

import os
import re
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext as _build_ext


def _numpy_include_dirs():
    # numpy is guaranteed available at build time via pyproject.toml
    import numpy as np
    inc = [np.get_include()]
    # keep your extra include path
    inc.append(str(Path(__file__).resolve().parent / "include"))
    return inc


class build_ext(_build_ext):
    """
    Build the C shared objects and then rename them to stable filenames:
      - strip ABI tags (e.g., .cpython-311-x86_64-linux-gnu)
      - force extension to .so (even on MacOS/Windows)
    """

    # Matches "...<name>.<abi>.<ext>" where ext is so/pyd/dll/dylib
    _abi_pat = re.compile(r"^(?P<stem>.+?)\.[^.]+\.(?P<ext>so|pyd|dll|dylib)$", re.IGNORECASE)

    def get_ext_filename(self, ext_name: str) -> str:
        # Let setuptools decide the default filename first (includes ABI tags usually)
        default = super().get_ext_filename(ext_name)
        # Then normalize it: strip ABI tag + force .so suffix
        base = os.path.basename(default)
        m = self._abi_pat.match(base)
        if m:
            base = f"{m.group('stem')}.so"
        else:
            # Fallback: replace final suffix with .so
            base = os.path.splitext(base)[0] + ".so"
        return os.path.join(os.path.dirname(default), base)

    def build_extension(self, ext: Extension) -> None:
        # Build to the *default* full path first
        super().build_extension(ext)

        # Figure out the file we just built (default name)
        built_path = Path(super().get_ext_fullpath(ext.name))
        # Figure out the normalized target path (.so, no ABI)
        target_path = Path(self.get_ext_fullpath(ext.name))

        if built_path == target_path:
            return

        target_path.parent.mkdir(parents=True, exist_ok=True)

        # If setuptools already placed a file at the target path, remove it
        if target_path.exists():
            target_path.unlink()

        built_path.replace(target_path)


ext_modules = [
    Extension("albopictus.modelAalbopictus03",
              ["src/albopictus/incubator03.c", "src/albopictus/modelAalbopictus03.c"]),
    Extension("albopictus.modelAalbopictus08",
              ["src/albopictus/gamma.c", "src/albopictus/incubator.c", "src/albopictus/modelAalbopictus08.c"]),
    Extension("albopictus.modelAalbopictus08b",
              ["src/albopictus/gamma.c", "src/albopictus/incubator.c", "src/albopictus/modelAalbopictus08b.c"]),
    Extension("albopictus.modelAalbopictus08c",
              ["src/albopictus/gamma.c", "src/albopictus/incubator.c", "src/albopictus/modelAalbopictus08c.c"]),
    Extension("albopictus.modelAalbopictus13",
              ["src/albopictus/gamma.c", "src/albopictus/incubator.c", "src/albopictus/modelAalbopictus13.c"]),
    Extension("albopictus.modelAalbopictus18",
              ["src/albopictus/ran_gen.c", "src/albopictus/gamma.c", "src/albopictus/spop.c",
               "src/albopictus/modelAalbopictus18.c"]),
    Extension("albopictus.modelCulex",
              ["src/albopictus/ran_gen.c", "src/albopictus/gamma.c", "src/albopictus/spop.c",
               "src/albopictus/modelCulex.c"]),
    Extension("albopictus.modelStochCHIKV",
              ["src/albopictus/ran_gen.c", "src/albopictus/spop01.c", "src/albopictus/gamma.c",
               "src/albopictus/modelStochCHIKV.c"]),
    Extension("albopictus.modelStochCDZ",
              ["src/albopictus/ran_gen.c", "src/albopictus/spop01.c", "src/albopictus/gamma.c",
               "src/albopictus/modelStochCDZ.c"]),
    Extension("albopictus.modelStochSand",
              ["src/albopictus/ran_gen.c", "src/albopictus/spop01.c", "src/albopictus/gamma.c",
               "src/albopictus/modelStochSand.c"]),
    Extension("albopictus.modelStochAalbopictus",
              ["src/albopictus/ran_gen.c", "src/albopictus/spop.c", "src/albopictus/gamma.c",
               "src/albopictus/modelStochAalbopictus.c"]),
    Extension("albopictus.modelStochAaegypti",
              ["src/albopictus/ran_gen.c", "src/albopictus/spop.c", "src/albopictus/gamma.c",
               "src/albopictus/modelStochAaegypti.c"]),
]

# apply include dirs to all extensions
for e in ext_modules:
    e.include_dirs = _numpy_include_dirs()

setup(
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
