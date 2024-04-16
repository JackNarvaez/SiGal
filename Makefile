PRINTMSG=--no-print-directory

.PHONY: setup

setup:
	@mkdir Data
	@cd setup;${MAKE} $(PRINTMSG) $@

galaxy:
	@cd src; ${MAKE} $(PRINTMSG) $@

animation3D:
	@cd plot; ${MAKE} $(PRINTMSG) $@

animation2D:
	@cd plot; ${MAKE} $(PRINTMSG) $@

clean:
	@cd src; ${MAKE} $(PRINTMSG) $@;\
	cd ../setup; ${MAKE} $(PRINTMSG) $@;\
	cd ../plot; ${MAKE} $(PRINTMSG) $@;
