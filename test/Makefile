.PHONY: all clean test
all: test

clean:
	@find .. -name '*~' -delete

test: ../ztables.py
	@./test.sh | grep ^test:

test-verbose: ../ztables.py
	@./test.sh
