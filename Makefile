# Builds and tests both solver apps; each app also has its own Makefile.
APPS := freesurf-2d freesurf-3d

.PHONY: all test debug clean run2d run3d

all test debug clean:
	@for app in $(APPS); do $(MAKE) -C $$app $@ || exit 1; done

run2d:
	$(MAKE) -C freesurf-2d run

run3d:
	$(MAKE) -C freesurf-3d run
