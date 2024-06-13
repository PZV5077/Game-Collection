/* Enable Audio */
let enableAudio = document.getElementById("enableAudio");
let body = document.getElementsByTagName("body")[0];

body.onclick = function () {
  body.removeChild(enableAudio);
  body.addEventListener("click", (e) => {
    e.stopPropagation();
  }, true);
};
