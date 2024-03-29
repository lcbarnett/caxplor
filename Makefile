WITH_GD       = 1
WITH_X11      = 1
WITH_PTHREADS = 1

SRC = main.c word.c ca.c rtab.c analyse.c sim_ana.c sim_bmark.c sim_test.c utils.c clap.c mt64.c strman.c

OBJ = $(patsubst %.c,.%.o,$(SRC))
DEP = $(patsubst %.o,%.d,$(OBJ))
BIN = caxplor

# Note: _GNU_SOURCE need for sincos function

OFLAGS  = -march=native -O3 -flto
WFLAGS  = -Wall -Werror -Wextra -Wconversion -Winline -Wno-unused-parameter
DFLAGS  = -D_DEFAULT_SOURCE -DUNSAFE_ZPIXMAP -D_GNU_SOURCE

CFLAGS  = $(OFLAGS) $(WFLAGS) $(DFLAGS)
LDFLAGS = $(OFLAGS) -lm

ifeq ($(WITH_GD),1)
	SRC     += cagd.c
	DFLAGS  += -DHAVE_GD
	LDFLAGS += -lgd
endif

ifeq ($(WITH_X11),1)
	SRC     += caX11.c sim_xplor.c screen_metrics.c
	DFLAGS  += -DHAVE_X11
	LDFLAGS += -lX11
endif

ifeq ($(WITH_PTHREADS),1)
	CC      += -pthread
	SRC     += sim_ddf.c sim_ddr.c
	DFLAGS  += -DHAVE_PTHREADS
	LDFLAGS += -lpthread
endif

REPDEP = sed -i -e '1s,\($*\)\.o[ :]*,\1.o \.$*.d: ,' \.$*.d

.PHONY: all clean diag

all: $(BIN)

clean:
	rm -f $(OBJ) $(DEP) $(BIN)

diag:
	@echo "*** SRC     = " $(SRC)
	@echo "*** OBJ     = " $(OBJ)
	@echo "*** DEP     = " $(DEP)
	@echo "*** BIN     = " $(BIN)
	@echo "*** WITH_GD = " $(WITH_GD)

$(OBJ): .%.o: %.c
	$(CC) -std=c99 -c -MMD -MP $(CFLAGS) $< -o $@
	@$(REPDEP)

$(BIN): $(OBJ)
	$(CC) $(OBJ) $(LDFLAGS) -o $(BIN)

-include $(DEP)
