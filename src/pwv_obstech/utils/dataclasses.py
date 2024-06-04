class Autocast:

    def __post_init__(self) -> None:

       for name, field in self.__dataclass_fields__.items():

            value = getattr(self, name)
            type = field.type

            if not isinstance(value, type):
                object.__setattr__(self, name, type(value))                      
    def __setattr__(self, name: str, value: object) -> None:

        fields = self.__dataclass_fields__
        if name in fields:
            value = fields[name].type(value)
