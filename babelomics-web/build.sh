#!/bin/sh
mkdir -p build
rm -rf build/index.html build/index.js build/workers build/fonts

vulcanize babelomics-index.html -o build/index.html --inline --strip --csp

cp -r bower_components/fontawesome/fonts build/
cp -r src/fonts/*.woff* build/fonts/

mkdir -p build/workers
cp lib/jsorolla/src/lib/worker-fileupload.js build/workers/

sed -i s@WORKERS_PATH=\"lib/jsorolla/src/lib/\"@WORKERS_PATH=\"workers/\"@g build/index.js

sed -i s@../bower_components/fontawesome/fonts/fontawesome-webfont.@fonts/fontawesome-webfont.@g build/index.html
sed -i s@../src/fonts/@fonts/@g build/index.html
