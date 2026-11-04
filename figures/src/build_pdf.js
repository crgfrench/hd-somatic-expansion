const {chromium} = require('playwright');
const fs = require('fs');
const TMP = process.argv[2];
(async () => {
  const pages = JSON.parse(fs.readFileSync(TMP + '/pages.json', 'utf8'));
  const b = await chromium.launch();
  for (let k = 0; k < pages.length; k++) {
    const [path, h] = pages[k];
    const p = await b.newPage({viewport: {width: 1520, height: h}});
    await p.goto('file://' + path);
    await p.waitForTimeout(2000);
    await p.pdf({path: `${TMP}/page${k}.pdf`, width: '1520px', height: h + 'px',
                 printBackground: true, pageRanges: '1',
                 margin: {top: '0', bottom: '0', left: '0', right: '0'}});
    await p.close();
  }
  await b.close();
})();
