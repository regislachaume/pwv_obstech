class Autocast:

    def __post_init__(self) -> None:

       for name, field in self.__dataclass_fields__.items():

            try: 
                value = getattr(self, name)
            except AttributeError:
                continue

            type_ = field.type

            if not isinstance(value, type_):
                object.__setattr__(self, name, type_(value))                      

    def __setattr__(self, name: str, value: object) -> None:

        fields = self.__dataclass_fields__
        
        if name in fields:

            type_ = fields[name].type

            if not isinstance(value, type_):
                value = type_(value)

        object.__setattr__(self, name, value)
