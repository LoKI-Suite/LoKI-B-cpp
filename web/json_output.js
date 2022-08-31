/* This function is called at the end of run. The arguments are a C-style string
 * (a character array) and its length. It represents the results of calling dump()
 * on the JSON object that is used in the C++ code to hold the output data.
 */
export function handleJsonOutput(str /*msg_ptr, msg_size*/) {
  // const msg = new Uint8Array(module.HEAPU8.buffer, msg_ptr, msg_size);
  // const str = String.fromCharCode.apply(null, msg);
  // const str = msg.reduce(
  //   (tot, i) => (tot += String.fromCharCode.apply(null, [i])),
  //   ""
  // );
  const json = JSON.parse(str);
  console.log(json);
  document.getElementById("json-results").innerHTML = JSON.stringify(
    json,
    {},
    2
  );
}
