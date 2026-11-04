// Evaluate Main.dc.html's renderVals() outside the canvas runtime and print its
// result, so the standalone page can carry the disease table as plain markup.
const fs = require('fs');
const path = require('path');
const P = process.argv[2] || path.join(__dirname, 'Main.dc.html');
const s = fs.readFileSync(P, 'utf8');

const open = s.indexOf('renderVals()');
if (open < 0) throw new Error('renderVals not found');
let i = s.indexOf('{', open), depth = 0, end = -1;
for (let k = i; k < s.length; k++) {
  if (s[k] === '{') depth++;
  else if (s[k] === '}') { depth--; if (depth === 0) { end = k; break; } }
}
const body = s.slice(i + 1, end);

// the component's props, as declared in data-props, with their defaults
const m = s.match(/data-props='([^']*)'/);
const props = {};
if (m) for (const [k, v] of Object.entries(JSON.parse(m[1]))) {
  if (k[0] !== '$' && v && 'default' in v) props[k] = v.default;
}

const self = { props };
const fn = new Function('return function(){' + body + '}')();
process.stdout.write(JSON.stringify(fn.call(self)));
