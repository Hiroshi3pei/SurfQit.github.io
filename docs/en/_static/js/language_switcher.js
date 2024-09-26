function getCurrentLanguageAndPath() {
    const path = window.location.pathname;
    const parts = path.split('/');
    
    // 'docs' を探してそのインデックスを取得
    const docsIndex = parts.indexOf('SurfQit.github.io'); //'docs' 'SurfQit.github.io'
    
    // 英語版かどうかを判定
    const isEnglish = parts[docsIndex + 1] === 'en';

    // basePath を設定
    const basePath = parts.slice(0, docsIndex + 1).join('/');

    // relativePath を設定
    let relativePath;
    if (isEnglish) {
        // 'docs/en/' の後の部分を取得
        relativePath = parts.slice(docsIndex + 2).join('/');
    } else {
        // 'docs/' の後の部分を取得
        relativePath = parts.slice(docsIndex + 1).join('/');
    }
    
    return {
        lang: isEnglish ? 'en' : 'ja',
        basePath: basePath,
        relativePath: relativePath || 'index.html'
    };
}

function setLanguageLinks() {
    const { lang, basePath, relativePath } = getCurrentLanguageAndPath();

    // 日本語と英語のリンクを動的に設定
    const japaneseLink = document.getElementById('japanese-link');
    const englishLink = document.getElementById('english-link');

    if (lang === 'en') {
        // 英語版から日本語版へ: 'en' を削除したパスに変更
        japaneseLink.href = `${basePath}/${relativePath}`;
        englishLink.href = `${basePath}/en/${relativePath}`;
    } else {
        // 日本語版から英語版へ: 'en' を追加したパスに変更
        japaneseLink.href = `${basePath}/${relativePath}`;
        englishLink.href = `${basePath}/en/${relativePath}`;
    }
}

document.addEventListener('DOMContentLoaded', function() {
    setLanguageLinks();
});