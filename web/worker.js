importScripts("loki_bindings.js");

async function lxcat_get(module, input) {
  for (const [i, file] of input.electronKinetics.LXCatFiles.entries()) {
    if (file.substr(0, 5) == "http:") {
      const response = await fetch(
        file, 
        {
          cache: "no-cache", 
          credentials: "include", 
          headers: {
            "Content-Type": "application/json", 
            "Authorization": "Bearer <TOKEN>"
          },
        }
      )

      // const valid = validate(JSON.parse(http.responseText));
      // if (!valid) console.log("AJV error: ", validate.errors);

      input.electronKinetics.LXCatFiles[i] = "/" + i + ".json";
      module.FS.writeFile("/" + i + ".json", await response.text(), {
        flags: "w+",
      });
    }
  }
}

async function worker_plot(x_ptr, x_size, y_ptr, y_size) {
  const x_arr = new Float64Array(this.HEAPF64.buffer, x_ptr, x_size);
  const y_arr = new Float64Array(this.HEAPF64.buffer, y_ptr, y_size);
  const data = new Array(x_size).fill().map((_, index) => {
    return { x: x_arr[index], y: y_arr[index] };
  });

  postMessage({ type: "DATA", results: data });
}

async function json_output(str) {
  postMessage({ type: "DONE", results: str });
}

self.onmessage = ({ data }) => {
  if (data.type === "COMPUTE") {
    create_module().then(async (module) => {
      parsed = JSON.parse(data.input);
      await lxcat_get(module, parsed);
      module.run(JSON.stringify(parsed), worker_plot.bind(module), json_output);
    });
  }
};
