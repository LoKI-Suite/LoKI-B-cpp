/* This function is called by the C++ function handleResults
 * (see loki-b/web/bindings.cpp). The arguments are the addresses
 * and sizes of two double-valued arrays that represent the energies
 * of LoKI-B's energy grid and the corresponsing values of the
 * EEDF in these points.
 *
 * NOTE: both arrays must have equal lengths.
 */
// export function plot(x_ptr, x_size, y_ptr, y_size) {
export function plot(data) {
  // Map sections of the emscripten heap onto contiguous arrays of doubles.
  // const x_arr = new Float64Array(module.HEAPF64.buffer, x_ptr, x_size);
  // const y_arr = new Float64Array(module.HEAPF64.buffer, y_ptr, y_size);
  // const data = new Array(x_size).fill().map((_, index) => {
  //   return { x: x_arr[index], y: y_arr[index] };
  // });
  const spec = {
    $schema: "https://vega.github.io/schema/vega-lite/v4.json",
    width: 800,
    height: 600,
    data: {
      values: data,
    },
    encoding: {
      x: { field: "x", type: "quantitative", title: "Energy (eV)" },
      y: {
        field: "y",
        type: "quantitative",
        scale: { type: "log" },
        title: "Eedf (eV^{-3/2})",
      },
    },
    mark: {
      type: "circle",
      // Enable tooltip so on mouseover it shows all data of that iteration
      tooltip: { content: "data" },
    },
    // Enable zooming and panning
    selection: { grid: { type: "interval", bind: "scales" } },
  };
  vegaEmbed(document.getElementById("plotting-area"), spec);
}
