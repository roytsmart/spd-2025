import pathlib
import spd2025


def test_all():

    paths = spd2025.figures.all()

    for path in paths:
        assert isinstance(path, pathlib.Path)
