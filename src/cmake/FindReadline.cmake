FIND_PATH(READLINE_INCLUDE_DIRS readline/readline.h)
FIND_LIBRARY(READLINE_LIBRARIES NAMES readline) 

FIND_PACKAGE_HANDLE_STANDARD_ARGS(READLINE DEFAULT_MSG READLINE_LIBRARIES READLINE_INCLUDE_DIRS)

