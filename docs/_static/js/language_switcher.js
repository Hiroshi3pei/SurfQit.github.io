// 現在のURLからベースパスを取得
function getBasePath() {
    const path = window.location.pathname;
    const parts = path.split('/');
    const langIndex = parts.indexOf('en');
    if (langIndex !== -1) {
        // 英語ページの場合
        return parts.slice(0, langIndex).join('/') + '/';
    } else {
        // 日本語ページの場合
        return parts.slice(0, parts.length - 1).join('/') + '/';
    }
}

// 言語切り替えリンクを更新
function updateLanguageLinks() {
    const basePath = getBasePath();
    const currentPath = window.location.pathname;
    const isEnglish = currentPath.includes('/en/');

    const japaneseLink = document.getElementById('japanese-link');
    const englishLink = document.getElementById('english-link');

    if (isEnglish) {
        japaneseLink.href = currentPath.replace('/en/', '/');
        englishLink.href = currentPath;
    } else {
        japaneseLink.href = currentPath;
        englishLink.href = basePath + 'en' + currentPath.substring(basePath.length - 1);
    }
}

// ページ読み込み時に実行
document.addEventListener('DOMContentLoaded', updateLanguageLinks);