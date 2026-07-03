"""Handle different kinds of error. To be expanded in future releases"""

class FileTypeError(TypeError):
    pass

class FileDataError(Exception):
    pass

class DirNotFoundError(Exception):
    pass
