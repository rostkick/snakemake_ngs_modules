class UserInput:
    __match_args__ = ('argument_name', 'argument_value')

    def __init__(self, argument_name: str, argument_value: str | bool | dict) -> None:
        self.argument_name = argument_name
        self.argument_value = argument_value


class CommandElement:

	def __init__(self, user_input: UserInput) -> None:
		self.user_input = user_input
		self.elements = []
		self.__indent_size = 4
		self.__sep = ' '
		self.__sep1 = '='
		self.__l_end = ' \\'

	@property
	def indent_size(self):
		return self.__indent_size

	@indent_size.setter
	def indent_size(self, indent_size) -> None:
		self.__indent_size = indent_size

	def __str(self, indent):
		lines = []
		ind = ' ' * (indent * self.__indent_size)
		ind1 = ' ' * ((indent + 1) * self.__indent_size)
		match self.user_input:
			case UserInput(str(argument_name), str(argument_value)):
				lines.append(f'{ind}{argument_name}{self.__sep}{argument_value}{self.__l_end}')
			case UserInput(str(argument_name), argument_value=True):
				lines.append(f'{ind}{argument_name}{self.__l_end}')
			case UserInput(str(argument_name), dict(argument_value)):
				lines.append(f'{ind}{argument_name}{self.__l_end}')
				for k, val in argument_value.items():
					lines.append(f'{ind1}{k}{self.__sep1}{val}{self.__l_end}')
			case _:
				pass

		for e in self.elements:
			lines.append(e.__str(indent + 1))

		lines = list(filter(None, lines))

		return '\n'.join(lines)

	def __str__(self):
		return self.__str(0).removesuffix(self.__l_end)


class CommandBuilder:
	def __init__(self, root_name) -> None:
		self.root_name = root_name
		self.__root = CommandElement(root_name)

	def add_cmd_element(self, user_input):
		self.__root.elements.append(
			CommandElement(user_input)
		)

	def set_indent_size(self, indent_size) -> None:
		for e in self.__root.elements:
			e.indent_size = indent_size

	def __str__(self) -> str:
		return str(self.__root)

class SlurmDecorator:
	def __init__(self, wrapped_builder: CommandBuilder, slurm_builder: CommandBuilder) -> None:
		self.slurm_builder = slurm_builder
		self.wrapped_builder = wrapped_builder

	def __str__(self) -> str:
		return '\n\n'.join((str(self.slurm_builder),
						str(self.wrapped_builder),
						'ENDINPUT'))
