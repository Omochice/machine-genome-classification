from os import PathLike


class MyData:
    def __init__(self, source_csv: PathLike) -> None:
        self._raw_data = ""
        self.params = {}
        self.images = []

        pass

    def _preprocess():
        pass

    def dump(self, dst: PathLike) -> None:
        pass

    def split()
