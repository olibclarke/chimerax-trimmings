"""The main startup script owns shapealign without a separate helper file."""

import importlib.util
from pathlib import Path
import runpy
import sys
import unittest
from unittest.mock import patch

from test_shape_align import cloud, model, session


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "chimerax_trimmings.py"


class MainScriptIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Initialize native settings/parsers before mocking registration.
        cls.s = session()

    def setUp(self):
        self.s.models.close(self.s.models.list())
        self.s.undo.clear()

    def test_import_is_safe_without_a_session_or_sibling_helper(self):
        spec = importlib.util.spec_from_file_location("trimmings_import_test", SCRIPT)
        module = importlib.util.module_from_spec(spec)
        with patch.dict(sys.modules, {"shape_align": None}), \
                patch("chimerax.core.commands.register") as register:
            spec.loader.exec_module(module)
        register.assert_not_called()
        self.assertTrue(callable(module.shape_align))
        self.assertTrue(callable(module.next_model))
        self.assertTrue(callable(module.register_command))
        self.assertFalse((ROOT / "shape_align.py").exists())

    def test_loading_main_script_registers_new_and_existing_commands(self):
        s = self.s
        with patch.dict(sys.modules, {"shape_align": None}), \
                patch("chimerax.core.commands.register") as register:
            runpy.run_path(str(SCRIPT), init_globals={"session": s})
        names = [call.args[0] for call in register.call_args_list]
        for name in ("shapealign", "nextmodel", "prevmodel", "hydrophobicpatches",
                     "modelmotion active", "reorient", "cyclemodeldisplay"):
            self.assertIn(name, names)
        self.assertEqual(names.count("shapealign"), 1)

    def test_native_shapealign_works_after_loading_only_main_script(self):
        from chimerax.core.commands import run

        s = self.s
        source = model(s, "source", cloud()[::5])
        target = model(s, "target", cloud()[::5])
        with patch.dict(sys.modules, {"shape_align": None}):
            runpy.run_path(str(SCRIPT), init_globals={"session": s})
            _, stats = run(s, f"shapealign #{source.id_string} to #{target.id_string} "
                          "objective gradient fallback_correlation 0")
        self.assertEqual(stats["objective"], "gradient")
        self.assertFalse(stats["partial"])
        self.assertGreater(stats["score"], .95)

    def test_startup_has_only_the_two_required_scripts(self):
        lines = (ROOT / "chimerax_trimmings.txt").read_text().splitlines()
        loaders = [line for line in lines if line.startswith("runscript ")]
        self.assertEqual(len(loaders), 2)
        self.assertTrue(loaders[0].endswith('chimerax_trimmings.py"'))
        self.assertTrue(loaders[1].endswith('map_salami.py"'))
        self.assertNotIn("shape_align.py", "\n".join(loaders))


if __name__ == "__main__":
    unittest.main()
