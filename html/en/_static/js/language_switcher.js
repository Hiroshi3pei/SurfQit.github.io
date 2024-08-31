// _static/language_switcher.js
document.addEventListener("DOMContentLoaded", function() {
  var selector = document.getElementById("language-selector").getElementsByTagName("select")[0];
  selector.addEventListener("change", function() {
    var lang = this.value;
    var currentUrl = window.location.pathname;
    var newUrl = currentUrl.replace(/\/(ja|en)\//, lang + "/");
    window.location.href = newUrl;
  });
});