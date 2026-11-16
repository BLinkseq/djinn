VERSION=3.0

djinn:
	@echo "Building djinn"
	mkdir -p bin/
	CGO_ENABLED=0 go build -ldflags "-X github.com/blinkseq/djinn/cmd.VERSION=$(VERSION) -s -w" -o bin/$@
