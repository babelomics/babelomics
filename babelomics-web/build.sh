#!/bin/sh
mkdir -p build
rm -rf build/index.html build/index.js build/fonts build/images

#vulcanize babelomics-index.html -o build/index.html --strip --csp --inline --config vulcanize.json
vulcanize babelomics-index.html -o build/index.html --inline --strip --csp --config vulcanize.json

cp -r bower_components/fontawesome/fonts build/
cp -r src/fonts/*.woff* build/fonts/
cp -r src/images build/
cp -r lib/jsorolla/styles/img/* build/images/
cp babelomics-config.js build/

sed -i s@../bower_components/fontawesome/fonts/fontawesome-webfont.@fonts/fontawesome-webfont.@g build/index.html
sed -i s@../src/fonts/@fonts/@g build/index.html
sed -i s@../src/images/@images/@g build/index.html

sed -i s@../lib/jsorolla/styles/img/@images/@g build/index.html
sed -i s@../lib/jsorolla/styles/fonts/@fonts/@g build/index.html

sed -i s@../babelomics-config.js@babelomics-config.js@g build/index.html
