#EnergyKey **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    Makefile                                           :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: uahmed <uahmed@student.hive.fi>            +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2024/07/14 20:03:19 by uahmed            #+#    #+#              #
#    Updated: 2024/07/14 20:03:27 by uahmed           ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

NAME        :=  scf
ERRTXT      :=  error.txt
OBJSDIR     :=  build
INCSDIR	    :=	include
SRCSDIR	    :=	src
DEPSDIR     :=  .deps

RM          :=  rm -rf
AR          :=  ar -rcs
CC          :=  c++
CFLAGS      :=  -Wall -Werror -Wextra -std=c++17
DEBUGFLAGS  :=  -g -fsanitize=address
DEPFLAGS    :=  -c -MT $$@ -MMD -MP -MF $(DEPSDIR)/$$*.d
SCREENCLR   :=  printf "\033c"
SLEEP       :=  sleep .1

SOURCES     :=  main.cpp \
		HF_SCF.cpp \
		MP2.cpp \
		helpers.cpp \

SRCS        :=  $(foreach file, $(SOURCES), $(shell find $(SRCSDIR) -name $(file)))
OBJS        :=  $(patsubst $(SRCSDIR)/%.cpp, $(OBJSDIR)/%.o, $(SRCS))
DEPS        :=  $(patsubst $(SRCSDIR)/%.cpp, $(DEPSDIR)/%.d, $(SRCS))
INCS	    :=	$(foreach header, $(INCSDIR), -I $(header))
INCS	    :=	$(RAYLIBINCS)

F           =   =====================================
B           =   \033[1m
T           =   \033[0m
G           =   \033[32m
V           =   \033[35m
C           =   \033[36m
R           =   \033[31m
Y           =   \033[33m

vpath %.cpp $(SRCSDIR)

define cc_cmd
$(OBJSDIR)/%.o: %.cpp | $(OBJSDIR) $(DEPSDIR)
	@if ! $(CC) -g $(INCS) $(DEPFLAGS) $$< -o $$@ 2> $(ERRTXT); then \
		printf "$(R)$(B)\nERROR!\n$(F)$(T)\n"; \
		printf "$(V)Unable to create object file:$(T)\n\n"; \
		printf "$(R)$(B)$$@$(T)\n"; \
		printf "$(Y)\n"; sed '$$d' $(ERRTXT); \
		printf "$(R)$(B)\n$(F)\nExiting...$(T)\n"; exit 1 ; \
	else \
		printf "$(C)$(B)☑$(T)$(V) $$<$ \n    $(C)⮑\t$(G)$(B)$$@$(T) \t\n"; \
	fi
endef

all: title $(NAME)

$(NAME): $(OBJS)
	@$(CC) -g $(INCS) -o $@ $^
	@make finish

debug: CFLAGS += $(DEBUGFLAGS)
debug: all

clean:
	@$(RM) $(OBJSDIR) $(DEPSDIR) $(ERRTXT)

fclean: clean
	@$(RM) $(NAME)

re: fclean all

title:
	@$(SCREENCLR) && printf "\n"
	@printf "$(C)╔╦╗╔═╗╔╦╗╦ ╦╦  ╔═╗╔═╗╔═╗$(T)\n"
	@printf "$(C)║║╠║ ║ ║║║ ║║  ║╣ ║║║║║║$(T)\n"
	@printf "$(C)╩ ╩╚═╝═╩╝╚═╝╩═╝╚═╝╚═╝╚═╝  by uahmed$(T)\n"
	@printf "$(G)$(B)$(F)\n$(T)\n"

finish:
	@printf "\n$(G)$(B)$(F)$(T)\n"
	@printf "$(C)╔═╗╦╔╗╔╦╔═╗╦ ╦╔═╗╔╦╗        $(V)$(B)$(NAME)$(T)\n"
	@printf "$(C)╠╣ ║║║║║╚═╗╠═╣║╣  ║║$(T)\n"
	@printf "$(C)╚  ╩╝╚╝╩╚═╝╩ ╩╚═╝═╩╝$(T)\n\n"

$(OBJSDIR) $(DEPSDIR):
	@mkdir -p $@

-include $(DEPS)

$(foreach build, $(OBJSDIR), $(eval $(call cc_cmd, $(build))))

.PHONY: all debug clean fclean re title finish
