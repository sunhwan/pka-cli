import os.path as osp

import pytest

from pka_cli import _env
from pka_cli.core import PkaPredictor


def _molgpka_available():
    return osp.isdir(osp.join(_env.molgpka_home(), "src"))


@pytest.fixture(scope="session")
def predictor():
    if not _molgpka_available():
        pytest.skip(
            "MolGpKa not found; run scripts/setup_molgpka.sh or set MOLGPKA_HOME"
        )
    return PkaPredictor()
