-- Let inline code wrap at punctuation in the LaTeX/PDF build.
local specials = {
  ['\\'] = '\\textbackslash{}', ['{'] = '\\{', ['}'] = '\\}', ['_'] = '\\_',
  ['$'] = '\\$', ['&'] = '\\&', ['#'] = '\\#', ['%'] = '\\%',
  ['^'] = '\\textasciicircum{}', ['~'] = '\\textasciitilde{}',
}
local breakable = { ['_'] = true, ['/'] = true, [','] = true, ['='] = true,
  ['.'] = true, ['-'] = true, [':'] = true, ['|'] = true }

local function wrap(text)
  local out = {}
  local prev = ''
  for ch in text:gmatch('.') do
    if prev:match('%l') and ch:match('%u') then out[#out + 1] = '\\allowbreak{}' end
    prev = ch
    out[#out + 1] = specials[ch] or ch
    if breakable[ch] then out[#out + 1] = '\\allowbreak{}' end
  end
  return '\\texttt{' .. table.concat(out) .. '}'
end

function Header(h) return h, false end
function Code(c)
  if FORMAT:match('latex') or FORMAT:match('pdf') then
    return pandoc.RawInline('latex', wrap(c.text))
  end
end
