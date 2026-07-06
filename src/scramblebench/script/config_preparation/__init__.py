"""This entails all of the key and value of the user's YAML as separate class. The reason why
default class is used instead of dataclass is that I want to have more control of the data that
the user used. Mainly for validating the keys and values. As each config has their own validation methods,
I could not use dependency injection like `validate(ConfigClass). Hence, I have decided to use the other way
(ConfigClass.validate_config())."""