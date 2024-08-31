// _static/language_switcher.js
document.addEventListener("DOMContentLoaded", function() {
  var selector = document.getElementById("language-selector").getElementsByTagName("select")[0];
  selector.addEventListener("change", function() {
    var lang = this.value;
    var currentUrl = window.location.pathname;
    
    // ルートにいる場合のパス設定
    var newUrl;
    if (currentUrl.endsWith("index.html") || currentUrl === "/") {
      if (lang === "en") {
        newUrl = "./en/index.html";
      } else {
        newUrl = "./index.html"; // 日本語はルートディレクトリ
      }
    } else {
      // jaまたはenのフォルダ内にいる場合のパス設定
      newUrl = currentUrl.replace(/\/(ja|en)\//, "/" + lang + "/");
    }
    
    window.location.href = newUrl;
  });
});