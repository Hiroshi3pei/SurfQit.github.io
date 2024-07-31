function initCopyButtons() {
  const buttons = document.querySelectorAll('.copy-button');
  buttons.forEach(button => {
    button.addEventListener('click', function() {
      const targetId = this.getAttribute('data-target');
      copyButton(targetId);
    });
  });
}

function copyButton(elementId) {
  var element = document.getElementById(elementId);
  navigator.clipboard.writeText(element.textContent);
}

// DOMの読み込みが完了したら初期化関数を呼び出す
document.addEventListener('DOMContentLoaded', initCopyButtons);