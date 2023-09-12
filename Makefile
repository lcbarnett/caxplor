SRC = main.c word.c ca.c rtab.c screen_metrics.c sim_xplor.c sim_ana.c sim_bmark.c sim_test.c utils.c clap.c mt64.c strman.c
OBJ = $(patsubst %.c,.%.o,$(SRC))
DEP = $(patsubst %.o,%.d,$(OBJ))
BIN = caxplor

# Note: _GNU_SOURCE need for sincos function

OFLAGS = -march=native -O3 -flto
WFLAGS = -Wall -Werror -Wextra -Wconversion -Winline -Wno-unused-parameter
DFLAGS = -D_DEFAULT_SOURCE -DUNSAFE_ZPIXMAP -D_GNU_SOURCE -DDFT_SINGLE_PREC_FLOAT
# DFLAGS = -D_DEFAULT_SOURCE -DUNSAFE_ZPIXMAP -D_GNU_SOURCE
CFLAGS = $(OFLAGS) $(WFLAGS) $(DFLAGS)

LDFLAGS = $(OFLAGS) -lgd -lX11 -lm

REPDEP = sed -i -e '1s,\($*\)\.o[ :]*,\1.o \.$*.d: ,' \.$*.d

.PHONY: all clean diag

all: $(BIN)

clean:
	rm -f $(OBJ) $(DEP) $(BIN)

diag:
	@echo "*** SRC = " $(SRC)
	@echo "*** OBJ = " $(OBJ)
	@echo "*** DEP = " $(DEP)
	@echo "*** BIN = " $(BIN)

$(OBJ): .%.o: %.c
	$(CC) -std=c99 -c -MMD -MP $(CFLAGS) $< -o $@
	@$(REPDEP)

$(BIN): $(OBJ)
	$(CC) $(OBJ) $(LDFLAGS) -o $(BIN)

-include $(DEP)
