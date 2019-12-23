var Module = require('./loki');

Module.onRuntimeInitialized = function() {
    const fileName = 'default_lokib_setup.in';
    const buffer_size = Module.lengthBytesUTF8(fileName) + 1;

    let buffer = Module._malloc(buffer_size);

    const run = Module.cwrap('run', 'number', '[number]');

    Module.stringToUTF8(fileName, buffer, buffer_size);

    console.log(run([buffer]));

    Module._free(buffer);
};