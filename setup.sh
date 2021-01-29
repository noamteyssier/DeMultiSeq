#!/usr/bin/env bash

# installs cargo if not found
if ! hash cargo ; then
  curl https://sh.rustup.rs -sSf | sh
fi

# build bin directory
if [ ! -d "bin/" ]; then
  mkdir bin/
fi

# make optimized binary
cargo build --release

# symlink binary to bin
ln -s $(pwd)/target/release/demux bin/demux
