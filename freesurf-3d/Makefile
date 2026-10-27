# GNU GCC (Homebrew gcc-NN); override with `make CC=...`
ifeq ($(origin CC),default)
  CC := $(firstword $(foreach v,16 15 14 13,$(shell command -v gcc-$(v) 2>/dev/null)))
  ifeq ($(CC),)
    $(error GNU GCC not found: run `brew install gcc` or pass CC=...)
  endif
endif
OPENMP   ?= 1
CFLAGS   ?= -O2
DBGFLAGS ?= -O0 -g -fsanitize=address,undefined -fno-omit-frame-pointer
WARN     := -std=gnu11 -Wall -Wextra -Wpedantic -Wshadow
LDLIBS   := -lm
ifeq ($(OPENMP),1)
  OMPFLAGS := -fopenmp
endif

BUILD  := build
OBJDIR := $(BUILD)/obj
DBGDIR := $(BUILD)/debug

LIB_SRC := src/util.c src/config.c src/sim.c src/bc.c src/vof.c src/momentum.c src/poisson.c \
           src/init.c src/solver.c src/diagnostics.c src/vtk_io.c src/output.c
LIB_OBJ := $(LIB_SRC:src/%.c=$(OBJDIR)/%.o)
DBG_OBJ := $(LIB_SRC:src/%.c=$(DBGDIR)/%.o)
TESTS   := $(BUILD)/test_plic3d $(BUILD)/test_poisson3d $(BUILD)/test_vof3d

.PHONY: all debug test run clean

all: $(BUILD)/freesurf3d

$(BUILD)/freesurf3d: $(OBJDIR)/main.o $(LIB_OBJ)
	$(CC) $(CFLAGS) $(OMPFLAGS) -o $@ $^ $(LDLIBS)

$(OBJDIR)/%.o: src/%.c | $(OBJDIR)
	$(CC) $(WARN) $(CFLAGS) $(OMPFLAGS) -MMD -MP -c -o $@ $<

# debug build is serial: sanitizers and OpenMP runtimes do not mix well
debug: $(BUILD)/freesurf3d_debug

$(BUILD)/freesurf3d_debug: $(DBGDIR)/main.o $(DBG_OBJ)
	$(CC) $(DBGFLAGS) -o $@ $^ $(LDLIBS)

$(DBGDIR)/%.o: src/%.c | $(DBGDIR)
	$(CC) $(WARN) -Wno-unknown-pragmas $(DBGFLAGS) -MMD -MP -c -o $@ $<

$(BUILD)/test_%: tests/test_%.c $(LIB_OBJ)
	$(CC) $(WARN) $(CFLAGS) $(OMPFLAGS) -Isrc -o $@ $^ $(LDLIBS)

test: $(TESTS)
	@for t in $(TESTS); do echo "== $$t"; ./$$t || exit 1; done
	@echo "all tests passed"

run: $(BUILD)/freesurf3d
	./$(BUILD)/freesurf3d cases/standing_wave_3d.cfg

$(OBJDIR) $(DBGDIR):
	mkdir -p $@

clean:
	rm -rf $(BUILD)

-include $(LIB_OBJ:.o=.d) $(DBG_OBJ:.o=.d) $(OBJDIR)/main.d $(DBGDIR)/main.d
