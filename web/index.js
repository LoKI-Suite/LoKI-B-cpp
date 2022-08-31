import { plot } from "./plot.js";
import { handleJsonOutput } from "./json_output.js";

export let module;
let file = "";
let global_input = {};
let results = {};

const file_input = document.getElementById("file-input");
const run_button = document.getElementById("run-button");
const dl_button = document.getElementById("dl-button");
const ref_button = document.getElementById("ref-button");

const worker = new Worker("./worker.js");
// const validate = new Ajv().compile(schema);

file_input.addEventListener("change", (event) => {
  // const list = event.target.files;
  file = event.target.files[0];
});

worker.onmessage = ({ data }) => {
  if (data.type === "DATA") {
    plot(data.results);
  }
  if (data.type === "DONE") {
    // handleJsonOutput(data.results);
    results = data.results;
    dl_button.style.display = "inline-block";
    ref_button.style.display = "inline-block";
  }
};

run_button.addEventListener("click", (_event) => {
  // Hide 'Download references' button
  dl_button.style.display = "none";
  ref_button.style.display = "none";

  // Wait for input file to be loaded
  file
    .text()
    .then((input) => {
      global_input = JSON.parse(input);
      worker.postMessage({
        type: "COMPUTE",
        input,
      });
    })
    .catch((err) => console.log(err));
});

dl_button.addEventListener("click", (_event) => {
  download("results.json", results, "text/json");
});

ref_button.addEventListener("click", (_event) => {
  let references = "";
  for (const file of global_input.electronKinetics.LXCatFiles) {
    if (file.substr(0, 5) == "http:") {
      const http = new XMLHttpRequest();
      http.open("GET", file, false);
      http.send(null);
      const lxcat_input = JSON.parse(http.responseText);
      for (const id in lxcat_input.references) {
        references += lxcat_input.references[id];
      }
    }
    // const lxcat_input = JSON.parse(
    //   module.FS.readFile(file, { encoding: "utf8" })
    // );
  }
  download("references.bib", references);
});

function download(filename, object, type = "text/plain") {
  var a = document.createElement("a");
  var file = new Blob([object], { type });
  a.href = URL.createObjectURL(file);
  a.download = filename;
  a.click();
}
