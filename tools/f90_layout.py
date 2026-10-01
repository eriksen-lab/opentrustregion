#!/usr/bin/env python3
# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
Check, and optionally fix, the layout of Fortran statements.

Every statement is checked for the column it starts at, which follows from the nesting
of the blocks (modules, procedures, interfaces, derived types, do, if, select, ...) it
sits in, and for how it is broken into continuation lines and spaced. The issues
reported are:

- indent: the statement does not start at the column its block nesting gives
- align: a continuation line does not start where its anchor puts it
- overlong: a line is longer than the column limit
- spacing: operators, punctuation or keywords are not written as rules 7 and 9 say
- layout: the statement is not broken as the rules below say
- comment: a comment line does not start at the column rule 10 gives
- name: an end statement does not name its block as rule 9 says
- literal: a real literal has no kind or a d exponent (rule 9), which is not fixed
- whitespace, blank, gap, newline: a line ends in blanks, a blank line is superfluous
  or missing, or the file does not end with a newline (rule 11)

The layout of a statement is the best one under these rules:

1. No line is longer than the column limit.
2. The statement uses as few lines as possible, subject only to the first part of rule
   4, and among the layouts that do, breaks each line as late as possible, subject to
   the rest of rule 4.
3. A continuation line is indented relative to the innermost open anchor. Anchors are
   brackets, the assignment operators ``=`` and ``=>``, the ``::`` of a declaration,
   the ``only:`` of a use statement and the ``=`` or ``=>`` of a declaration
   initializer. Breaking directly after an anchor gives a hanging indent of four
   columns, relative to the line a bracket is on or to the start of the statement for
   the other anchors; breaking later aligns with the first thing after the anchor. A
   bracket whose first argument is a bracket group hanging on the same line hangs with
   it. The hanging indent of an initializer in a declaration of several entities is
   instead four columns relative to the names of the entities. Without an open anchor,
   continuation lines are indented by four columns.
4. Bracket groups are kept whole in three tiers. First, before the line count: an array
   constructor, an array section, or the shape of an array in a declaration or an
   allocate statement is not split if it fits whole on a line, and neither are the
   contents of a group holding a single argument or index hung after its opening
   bracket, nor a comparison (``==``, ``/=``, ``<``, ``<=``, ``>``, ``>=``, their dotted
   forms, ``.eqv.`` and ``.neqv.``) that fits whole on a continuation line with a
   hanging indent, which extends from the start of its left to the end of its right
   operand, as far as the nearest looser operator, comma, assignment or bracket on
   either side. Second, among the layouts with the fewest lines: a group containing no
   other bracket group apart from grouping parentheses is not split if it fits whole,
   and no group that fits whole has its contents hung after its opening bracket. A
   group fits whole if it fits where it stands, where a break directly before it would
   put it, or on a continuation line with a hanging indent. Where no layout keeps a
   tier, the one breaking it the fewest times is used. Third, the layouts splitting the
   fewest operands are preferred, where an operand is split by a break after an
   operator or comma binding tighter than another one in the same bracket group, and
   counts once for every such group. So a line is broken at ``.or.`` rather than inside
   the comparison next to it, and at a comma rather than inside an argument. An
   operator or comma also binds tighter than any operator or comma outside its bracket
   group, and so does a break directly after an opening bracket, so that a call is kept
   whole where breaking outside it costs no line. The condition of a single-line ``if``
   is looser than the statement it controls, and the commas between the entities of a
   declaration are looser than an initializer.
5. A trailing ``then`` that does not fit takes a hanging indent of four columns like
   any other continuation line. Breaking before it is looser than any break in the
   condition.
6. An array constructor given to ``reshape`` is written one matrix row per line when
   every row fits, where a row holds as many items as the extent filled first: the
   first extent of the shape, or the one ``order`` names first. The shape is resolved
   from integer literals, from integer named constants and, for ``shape(a)``, from the
   extents of ``a``, as declared in the same file, as long as every declaration of the
   name agrees. When the shape cannot be resolved or the rows do not fit, the
   constructor follows the normal rules.
7. Binary ``+``, ``-``, ``*`` and ``/``, the comparisons (``==``, ``/=``, ``<``,
   ``<=``, ``>``, ``>=`` and their dotted forms) and ``.and.``, ``.or.``, ``.eqv.`` and
   ``.neqv.`` have a blank on either side, ``**``, ``//`` and ``%`` have none, and
   ``.not.`` is followed by one blank while a unary ``+`` or ``-`` is followed by none.
   ``=>``, the ``=`` of an assignment or initializer and the ``::`` of a declaration
   have a blank on either side, as has the ``::`` after the type of an allocate
   statement or an array constructor; the ``=`` of a keyword argument, like any other
   ``=`` inside brackets, has none. A ``:`` inside brackets has no blank on either
   side, and one outside them, after ``only`` or the name of a construct, has one blank
   after it. A comma has no blank before and one after it, an opening bracket none
   after it and a closing bracket none before it, and outside string literals no two
   blanks follow each other, so that code is never aligned by padding. The keyword of a
   statement taking a condition or selector (``if``, ``else if``, ``do while``, ``do
   concurrent``, ``select case``, ``select type``, ``select rank``, ``case``, ``type
   is``, ``class is``, ``rank``, ``where``, ``elsewhere``, ``forall``, ``associate``)
   is followed by one blank before its opening bracket, and every other name, such as
   that of a procedure, an array, a type, an attribute or another statement, by none.
   An operator named in ``operator(...)`` has no blanks.
8. The format of an input/output statement (the second positional or the ``fmt``
   argument of ``read`` and ``write``, and the first item of ``print``) is written in
   single quotes, and every other string literal in double quotes. String literals
   joined by ``//`` are merged into one, and a string literal is only
   split after a blank between two words, into ``"... "// &`` and ``"..."`` on the next
   line, where the word starting the next line contains a letter or digit, so that a
   word of punctuation only, such as ``|`` or ``=``, stays with the word before it. A
   literal containing two blanks in a row, such as a row of a table, is taken to be
   aligned by hand: it is neither merged nor split, and only moved as a whole. A
   literal that fits whole on a continuation line is kept whole like an
   innermost bracket group under rule 4; one that has to be split anyway fills lines
   like any other text, so it may start on the line before.
9. Keywords, the logical constants and the dotted operators are written in lower case.
   End statements, ``else if``, ``select case``, ``select type``, ``select rank`` and
   ``do while`` are written as two words separated by one blank, as in ``end if``. The
   end statement of a module, submodule, program, subroutine or function names it, as
   in ``end subroutine solver``, while that of a derived type or interface does not.
   The length and kind of a character type are written as ``len=`` and ``kind=``, as is
   the kind of an intrinsic conversion (``real``, ``int``, ``nint``, ``logical``,
   ``cmplx``, ``aint``, ``anint``, ``ceiling``, ``floor``, ``char``, ``achar``). A
   real literal has a digit on either side of its point, as in ``0.5_rp``, has no point
   when a whole mantissa takes an exponent, as in ``1e-14_rp``, and is written with a
   lower case ``e`` and an exponent without a plus sign or leading zeros. Every real
   literal carries a kind and none has a ``d`` exponent; since adding a kind changes
   the precision, a literal without one is only reported.
10. A comment line between two statements starts at the column of the statement after
    it, so that a comment before ``else``, ``case`` or ``end`` starts where that
    statement does. Comment lines inside a continued statement move with it. A comment
    has one blank between its ``!`` and its text, unless it starts with ``!!`` or is a
    directive such as ``!$omp``, and a comment after code is preceded by two blanks.
11. No line ends in blanks, no file starts or ends with a blank line, no two blank lines
    follow each other, and a file ends with one newline. A subroutine or function is
    followed by one blank line, unless it is the last one in an interface block, and
    the ``contains`` of a program unit or procedure has one blank line before and after
    it, while that of a derived type has none.
12. In a declaration of several entities, an entity with an initializer that fits whole
    on a continuation line with a hanging indent is kept whole like a comparison under
    the first tier of rule 4, while one that does not starts a new line, as does the
    entity after it, unless no layout does so.

Statements with comments, ``;`` or a continuation line starting with ``&`` are never
reflowed: they keep their line breaks and comments, and only their indentation,
continuation indentation and the spacing and keywords of rules 7 and 9 are fixed. Their
line length is checked but cannot be fixed.

By default only statements touched by the uncommitted changes relative to a git
reference are checked, so that the check can be run after every edit without reporting
the whole code base. Files that git does not track count as changed throughout. Use
--all to check whole files. Run from a subdirectory, only the files below it are
checked.
"""

import argparse
import re
import subprocess
import sys
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

LIMIT = 88
HANG = 4


# lines and statements


def code_part(line: str) -> str:
    """
    this function returns the line with string literals blanked and a trailing comment
    removed, so that column positions stay valid for the original line
    """
    out = []
    quote = None
    for ch in line:
        if quote:
            out.append(ch if ch == quote else " ")
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
            out.append(ch)
        elif ch == "!":
            break
        else:
            out.append(ch)
    return "".join(out).rstrip()


def comment_text(comment: str) -> str:
    """
    this function returns a comment starting with ! with one blank between the ! and its
    text, leaving directives such as !$omp and comments starting with !! as they are
    """
    return re.sub(r"^!(?![!$])\s*(?=\S)", "! ", comment.rstrip())


def is_code(line: str) -> bool:
    """
    this function returns whether a line carries code, as opposed to being blank, a
    comment or a preprocessor directive
    """
    stripped = line.lstrip()
    return bool(stripped) and not stripped.startswith(("!", "#"))


def indent_of(line: str) -> int:
    """
    this function returns the column the text of a line starts at
    """
    return len(line) - len(line.lstrip())


@dataclass
class Statement:
    """
    this class holds the lines a statement spans, together with those of them that
    carry code
    """

    first: int  # index of the first line
    last: int  # index of the last line, inclusive
    code: List[int] = field(default_factory=list)  # indices of the code lines


def statements(lines: List[str]) -> List[Statement]:
    """
    this function groups code lines into statements, where comment and blank lines
    between the lines of a continued statement belong to it
    """
    result = []
    current = None
    for i, line in enumerate(lines):
        if not is_code(line):
            continue
        if current is None:
            current = Statement(i, i)
        current.code.append(i)
        current.last = i
        if not code_part(line).endswith("&"):
            result.append(current)
            current = None
    if current is not None:
        result.append(current)
    return result


@dataclass
class Joined:
    """
    this class holds a statement joined onto one line, together with where it is
    currently broken and whether it may be reflowed
    """

    text: str  # the statement on one line
    code: str  # the same with string literals blanked
    breaks: List[int]  # positions the statement is currently broken at
    reflowable: bool


def join_statement(lines: List[str], stmt: Statement) -> Joined:
    """
    this function joins the lines of a statement into one, recording where it is
    currently broken
    """
    text, code, breaks = "", "", []
    reflowable = stmt.last - stmt.first + 1 == len(stmt.code)
    for k, i in enumerate(stmt.code):
        line = lines[i]
        blanked = code_part(line)
        if line[len(blanked) :].strip():
            reflowable = False
        start = indent_of(line)
        seg_code, seg_text = blanked[start:], line[start : len(blanked)]
        if k < len(stmt.code) - 1:
            seg_code = seg_code[:-1].rstrip()
            seg_text = seg_text[: len(seg_code)]
        if seg_code.startswith("&"):
            reflowable = False
        if k > 0:
            # no blank after an opening bracket, before a closing one, or after a
            # concatenation written without blanks
            tight = code.endswith(("(", "[")) or seg_code.startswith((")", "]"))
            tight |= code.endswith("//") and not code.endswith(" //")
            text, code = text + ("" if tight else " "), code + ("" if tight else " ")
        text, code = text + seg_text, code + seg_code
        if k < len(stmt.code) - 1:
            breaks.append(len(code))
    if ";" in code:
        reflowable = False
    return Joined(text, code, breaks, reflowable)


# block indentation


PROCEDURE = re.compile(
    r"((pure|impure|elemental|recursive|non_recursive|module)\s+)*"
    r"((integer|real|logical|complex|character|type|class)\s*(\([^)]*\))?\s+)?"
    r"((pure|impure|elemental|recursive|non_recursive|module)\s+)*"
    r"(subroutine|function)\s+\w+"
)


def block_role(code: str) -> Optional[str]:
    """
    this function returns whether a statement closes a block, continues one (else,
    case, contains, ...), opens one, or none of these
    """
    c = code.strip().lower()
    c = re.sub(r"^\w+\s*:\s*(?=(do|if|select|block|associate)\b)", "", c)
    if re.match(
        r"end(\s|$)|end(do|if|select|where|forall|associate|block|interface|type|"
        r"module|subroutine|function|program)\b",
        c,
    ):
        return "close"
    if re.match(
        r"(else|elseif|elsewhere|case|contains)\b|class\s+(is|default)\b|type\s+is\b"
        r"|rank\s*(\([^=]*\)|default)\s*$",
        c,
    ):
        return "middle"
    if re.match(r"module\s+procedure\b", c):
        return None
    if re.match(r"(module|submodule|program)\b", c):
        return "open"
    if re.match(r"(abstract\s+)?interface\b", c) or PROCEDURE.match(c):
        return "open"
    if re.match(r"type\s*(,[^:]*)?::", c) or re.match(r"type\s+\w+\s*$", c):
        return "open"
    if re.match(r"do(\s|$)|block\s*$|critical\s*$|associate\s*\(|enum\s*,", c):
        return "open"
    if re.match(r"select\s+(case|type|rank)\b", c):
        return "open"
    if re.match(r"if\s*\(", c) and re.search(r"\bthen$", c):
        return "open"
    if re.match(r"(where|forall)\s*\(", c):
        depth = 0
        for col, ch in enumerate(c):
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth == 0:
                    return "open" if not c[col + 1 :].strip() else None
    return None


def expected_indents(codes: List[str]) -> List[int]:
    """
    this function returns the column every statement, given by its joined code, is
    expected to start at
    """
    level, result = 0, []
    for code in codes:
        role = block_role(code)
        if role == "close":
            level = max(level - 1, 0)
            result.append(HANG * level)
        elif role == "middle":
            result.append(HANG * max(level - 1, 0))
        else:
            result.append(HANG * level)
            if role == "open":
                level += 1
    return result


def block_rules(
    lines: List[str], stmts: List[Statement], codes: List[str]
) -> Tuple[Dict[int, str], Set[int], Set[int]]:
    """
    this function returns what rules 9 and 11 ask of the blocks of a file: for every
    statement ending a program unit, procedure, derived type or interface, the end
    statement naming the program unit or procedure but not the derived type or
    interface, and the lines that must be followed by one blank line and those that
    must be followed by none, where a procedure is followed by one, unless an end
    interface follows it, and the contains of a program unit or procedure is surrounded
    by one, whereas that of a derived type is surrounded by none
    """
    ends, one, none = {}, set(), set()
    # the kind and name of every open block, or None for a block of neither
    stack: List[Tuple[Optional[str], Optional[str]]] = []

    def neighbour(i: int, step: int) -> int:
        while 0 <= i < len(lines) and not lines[i].strip():
            i += step
        return i

    for n, (stmt, code) in enumerate(zip(stmts, codes)):
        c = code.strip()
        role = block_role(code)
        if role == "open":
            procedure = re.search(r"\b(subroutine|function)\s+(\w+)", c, re.IGNORECASE)
            unit = re.match(r"(module|program)\s+(\w+)\s*$", c, re.IGNORECASE)
            unit = unit or re.match(
                r"(submodule)\s*\([^)]*\)\s*(\w+)\s*$", c, re.IGNORECASE
            )
            if PROCEDURE.match(c.lower()) and procedure:
                stack.append((procedure.group(1).lower(), procedure.group(2)))
            elif unit:
                stack.append((unit.group(1).lower(), unit.group(2)))
            elif re.match(r"(abstract\s+)?interface\b", c, re.IGNORECASE):
                stack.append(("interface", None))
            elif re.match(r"type\b", c, re.IGNORECASE):
                stack.append(("type", None))
            else:
                stack.append((None, None))
        elif role == "middle" and re.match(r"contains\b", c, re.IGNORECASE):
            around = none if stack and stack[-1][0] == "type" else one
            before = neighbour(stmt.first - 1, -1)
            if before >= 0:
                around.add(before)
            around.add(stmt.last)
        elif role == "close" and stack:
            kind, name = stack.pop()
            if kind is None:
                continue
            if re.fullmatch(rf"end(\s*{kind}(\s+\S.*)?)?", c, re.IGNORECASE):
                ends[n] = f"end {kind} {name}" if name else f"end {kind}"
            if kind in ("subroutine", "function"):
                after = neighbour(stmt.last + 1, 1)
                if after < len(lines):
                    ending = re.match(
                        r"end\s*interface\b", lines[after].strip(), re.IGNORECASE
                    )
                    (none if ending else one).add(stmt.last)
    return ends, one, none


# rewrites of a joined statement


def literals(text: str) -> List[Tuple[int, int]]:
    """
    this function returns the first and last position of every string literal in a
    statement
    """
    spans, i = [], 0
    while i < len(text):
        if text[i] in "'\"":
            quote, j = text[i], i + 1
            while j < len(text):
                if text[j] == quote:
                    if text[j + 1 : j + 2] != quote:
                        break
                    j += 1
                j += 1
            if j < len(text):
                spans.append((i, j))
            i = j
        i += 1
    return spans


def aligned_by_hand(text: str, first: int, last: int) -> bool:
    """
    this function returns whether the string literal from first to last is aligned by
    hand, such as a row of a table, which it is when it contains two blanks in a row
    """
    return "  " in text[first + 1 : last]


def merge_literals(joined: Joined) -> Joined:
    """
    this function returns the statement with every concatenation of two string literals
    written as one literal, unless one of them is aligned by hand
    """
    text, breaks = joined.text, list(joined.breaks)
    merged = True
    while merged:
        merged = False
        spans = literals(text)
        for (a_start, a_end), (b_start, b_end) in zip(spans, spans[1:]):
            same = text[a_end] == text[b_start]
            aligned = aligned_by_hand(text, a_start, a_end)
            aligned = aligned or aligned_by_hand(text, b_start, b_end)
            if (
                same
                and not aligned
                and re.fullmatch(r"\s*//\s*", text[a_end + 1 : b_start])
            ):
                cut = b_start + 1 - a_end
                breaks = [
                    b if b <= a_end else b - cut
                    for b in breaks
                    if not a_end < b <= b_start
                ]
                text = text[:a_end] + text[b_start + 1 :]
                merged = True
                break
    code = code_part(text)
    if len(code) != len(text):
        return joined
    return Joined(text, code, breaks, joined.reflowable)


def next_content(code: str, i: int) -> int:
    """
    this function returns the first position from i on that is not a blank
    """
    while i < len(code) and code[i] == " ":
        i += 1
    return i


def operand_before(code: str, i: int) -> bool:
    """
    this function returns whether the code before position i, ignoring blanks, ends an
    operand, which makes an operator at i binary
    """
    j = i - 1
    while j >= 0 and code[j] == " ":
        j -= 1
    if j < 0:
        return False
    ch = code[j]
    # the sign of an exponent such as 1e-3 is not an operator
    if ch in "eEdD" and j > 0 and code[j - 1] in "0123456789.":
        k = j - 1
        while k >= 0 and code[k] in "0123456789.":
            k -= 1
        if k < 0 or not (code[k].isalnum() or code[k] == "_"):
            return False
    # a dotted operator such as .gt. ends in a period but is no operand
    if ch == "." and re.search(
        r"\.(and|or|eqv|neqv|not|eq|ne|lt|le|gt|ge)\.$", code[: j + 1], re.IGNORECASE
    ):
        return False
    return ch.isalnum() or ch in "_)]\"'."


def space_operators(code: str) -> List[Tuple[int, int, str]]:
    """
    this function returns the replacements putting a blank on either side of every
    binary +, -, * and /, and none around **
    """
    if re.match(
        r"\s*(real|integer|logical|complex|character)\s*\*", code, re.IGNORECASE
    ):
        return []
    edits = []

    def around(i: int, width: int, replacement: str) -> None:
        start, end = i, i + width
        while start > 0 and code[start - 1] == " ":
            start -= 1
        while end < len(code) and code[end] == " ":
            end += 1
        edits.append((start, end, replacement))

    for i, ch in enumerate(code):
        if ch not in "+-*/":
            continue
        if code.startswith("**", i) and operand_before(code, i):
            around(i, 2, "**")
            continue
        if ch == "*" and "*" in (code[i - 1 : i], code[i + 1 : i + 2]):
            continue
        # the default format of print and read
        if ch == "*" and re.search(r"\b(print|read)\s*$", code[:i], re.IGNORECASE):
            continue
        if ch == "/" and (
            code[i + 1 : i + 2] in ("/", "=", ")") or code[i - 1 : i] in ("/", "(")
        ):
            continue
        if operand_before(code, i):
            around(i, 1, f" {ch} ")
    return edits


# keywords and the spacing around punctuation


# the keywords written in lower case; Fortran does not distinguish an identifier of the
# same name from the keyword, so lowering it changes nothing
KEYWORDS = frozenset("""
    abstract allocatable allocate associate asynchronous backspace bind block call case
    character class close codimension common complex contains contiguous continue
    critical cycle data deallocate default deferred dimension do double elemental else
    elsewhere end entry enum enumerator equivalence error exit extends external final
    flush forall format function generic go if implicit import impure in inout inquire
    integer intent interface intrinsic is kind len logical module namelist
    non_intrinsic non_overridable none nopass nullify only open operator optional out
    parameter pass pointer precision print private procedure program protected public
    pure rank read real recursive result return rewind save select sequence stop
    submodule subroutine target then to type use value volatile wait where while write
    """.split())

# the blocks closed by an end statement written as two words
END_BLOCKS = (
    "if|do|select|where|forall|associate|block|critical|enum|interface|type|module|"
    "submodule|subroutine|function|program"
)

# the statements whose keyword is followed by one blank before its opening bracket,
# whereas every other name is directly followed by its opening bracket
SPACED_STATEMENTS = (
    "else if|if|do while|do concurrent|select case|select type|select rank|case|"
    "type is|class is|rank|where|elsewhere|forall|associate"
)


def apply_edits(joined: Joined, edits: List[Tuple[int, int, str]]) -> Joined:
    """
    this function returns the statement with the given non-overlapping replacements of
    positions start to end applied, moving its breaks along
    """
    edits = sorted(e for e in edits if joined.code[e[0] : e[1]] != e[2])
    if not edits:
        return joined
    text, code = joined.text, joined.code
    for start, end, replacement in reversed(edits):
        text = text[:start] + replacement + text[end:]
        code = code[:start] + replacement + code[end:]

    def move(b: int) -> int:
        shift = 0
        for start, end, replacement in edits:
            # an insertion at a break goes onto the line after it
            if end < b or start < end == b:
                shift += len(replacement) - (end - start)
            elif start < b:
                return start + len(replacement) + shift
        return b + shift

    return Joined(text, code, [move(b) for b in joined.breaks], joined.reflowable)


def closing(code: str, o: int) -> int:
    """
    this function returns the position of the bracket closing the one at o, or the end of
    the statement when it is not closed
    """
    depth = 0
    for i in range(o, len(code)):
        if code[i] in "([":
            depth += 1
        elif code[i] in ")]":
            depth -= 1
            if depth == 0:
                return i
    return len(code)


def statement_starts(code: str) -> List[int]:
    """
    this function returns where the keyword of a statement starts, after the name of a
    construct, together with where the statement controlled by a single-line if or where
    starts
    """
    start = next_content(code, 0)
    m = re.match(r"\w+\s*:(?!:)\s*", code[start:])
    if m and re.match(
        r"(do|if|select|block|associate|critical|forall|where)\b",
        code[start + m.end() :],
        re.IGNORECASE,
    ):
        start += m.end()
    starts = [start]
    m = re.match(r"(if|where)\s*\(", code[start:], re.IGNORECASE)
    if m:
        after = next_content(code, closing(code, start + m.end() - 1) + 1)
        if after < len(code) and not re.match(r"then\s*$", code[after:], re.IGNORECASE):
            starts.append(after)
    return starts


def lower_keywords(code: str) -> List[Tuple[int, int, str]]:
    """
    this function returns the replacements writing keywords, logical constants and
    dotted operators in lower case
    """
    edits = []
    for m in re.finditer(r"\b[A-Za-z]\w*\b", code):
        word = m.group(0)
        if word.lower() in KEYWORDS and word != word.lower():
            edits.append((m.start(), m.end(), word.lower()))
    for m in re.finditer(
        r"\.(true|false|not|and|or|eqv|neqv|eq|ne|lt|le|gt|ge)\.", code, re.IGNORECASE
    ):
        edits.append((m.start(), m.end(), m.group(0).lower()))
    return edits


def keyword_words(code: str) -> List[Tuple[int, int, str]]:
    """
    this function returns the replacements writing end statements, else if, select case,
    select type and do while as two words separated by one blank
    """
    edits = []
    for start in statement_starts(code):
        for pattern in (
            rf"(end)\s*({END_BLOCKS})\b(?=\s*(\w+\s*)?$)",
            r"(else)\s*(if)\s*\(",
            r"(select)\s*(case|type|rank)\s*\(",
            r"(do)\s+(while)\s*\(",
        ):
            m = re.match(pattern, code[start:], re.IGNORECASE)
            if m:
                words = f"{m.group(1).lower()} {m.group(2).lower()}"
                edits.append((start + m.start(1), start + m.end(2), words))
                break
    return edits


def keyword_brackets(code: str) -> List[Tuple[int, int, str]]:
    """
    this function returns the replacements putting one blank between the keyword of a
    statement taking a condition or selector (if, do while, select, case, where, ...)
    and its opening bracket, and none between any other name and its opening bracket
    """
    edits, spaced = [], set()
    for start in statement_starts(code):
        m = re.match(rf"({SPACED_STATEMENTS})(\s*)\(", code[start:], re.IGNORECASE)
        if m:
            edits.append((start + m.start(2), start + m.end(2), " "))
            spaced.add(start + m.end(2))
    for m in re.finditer(r"\b([A-Za-z]\w*)(\s+)\(", code):
        if m.end(2) not in spaced:
            edits.append((m.start(2), m.end(2), ""))
    return edits


def arguments(code: str, o: int) -> List[Tuple[int, Optional[str]]]:
    """
    this function returns where every argument of the bracket group opened at o starts,
    together with its keyword, or None for a positional argument
    """
    result = []
    for start, end in pieces(code, o):
        first = next_content(code, start)
        if first < end:
            m = re.match(r"(\w+)\s*=(?![=>])", code[first:end])
            result.append((first, m.group(1).lower() if m else None))
    return result


# the intrinsic conversions whose second argument is the kind, and cmplx, whose third is
CONVERSIONS = "real|int|nint|logical|aint|anint|ceiling|floor|char|achar"


def keyword_arguments(code: str) -> List[Tuple[int, int, str]]:
    """
    this function returns the replacements writing out the len and kind of a character
    type and the kind of an intrinsic conversion
    """
    edits = []
    for m in re.finditer(
        rf"\b(character|cmplx|{CONVERSIONS})\s*\(", code, re.IGNORECASE
    ):
        name, args = m.group(1).lower(), arguments(code, m.end() - 1)
        if name == "character":
            names = {0: "len", 1: "kind"}
        else:
            names = {2: "kind"} if name == "cmplx" else {1: "kind"}
        for k, (first, keyword) in enumerate(args):
            if keyword is None and k in names:
                edits.append((first, first, f"{names[k]}="))
    return edits


# a real literal: a mantissa with a point or an exponent, and an optional kind, where a
# point directly followed by a dotted operator such as 1.eq.2 belongs to the operator
REAL_LITERAL = re.compile(
    r"(?<![\w.])(\d+\.(?!(eq|ne|lt|le|gt|ge|and|or|not|eqv|neqv|true|false)\.)\d*"
    r"|\.\d+|\d+(?=[eEdD][+-]?\d))([eEdD][+-]?\d+)?(_\w+)?(?![\w.])",
    re.IGNORECASE,
)


def real_literals(code: str) -> List[Tuple[int, int, str]]:
    """
    this function returns the replacements writing real literals with a digit on either
    side of a point, without a point in a whole mantissa with an exponent, and with a
    lower case e and an exponent without a plus sign or leading zeros; literals with a
    d exponent are left to the check for literals without a kind
    """
    edits = []
    for m in REAL_LITERAL.finditer(code):
        mantissa, exponent, kind = m.group(1), m.group(3), m.group(4) or ""
        if exponent and exponent[0] in "dD":
            continue
        whole, _, fraction = mantissa.partition(".")
        whole = whole or "0"
        if exponent:
            sign = "-" if "-" in exponent else ""
            digits = exponent[1:].lstrip("+-").lstrip("0") or "0"
            if fraction.strip("0"):
                wanted = f"{whole}.{fraction}e{sign}{digits}{kind}"
            else:
                wanted = f"{whole}e{sign}{digits}{kind}"
        else:
            wanted = f"{whole}.{fraction or '0'}{kind}"
        edits.append((m.start(), m.end(), wanted))
    return edits


def kindless_reals(code: str) -> List[str]:
    """
    this function returns the real literals without a kind or with a d exponent, whose
    precision only a hand edit can make explicit
    """
    return [
        m.group(0)
        for m in REAL_LITERAL.finditer(code)
        if not m.group(4) or (m.group(3) and m.group(3)[0] in "dD")
    ]


def unary_signs(code: str) -> List[Tuple[int, int, str]]:
    """
    this function returns the replacements removing the blanks after a unary + or -
    """
    edits = []
    for m in re.finditer(r"[+-](\s+)", code):
        if not operand_before(code, m.start()):
            edits.append((m.start(1), m.end(1), ""))
    return edits


def space_punctuation(code: str) -> List[Tuple[int, int, str]]:
    """
    this function returns the replacements putting one blank on either side of ::, =>,
    a comparison and the = of an assignment or initializer, and after a comma, none
    before a comma, around //, around the = of a keyword argument or any other = inside
    brackets, after an opening or before a closing bracket
    """
    edits, depth = [], 0
    pattern = (
        r"([(\[])\s*|\s*([)\]])|\s*(::|=>|//|==|/=|<=|>=|<|>|,|%|:|"
        r"(?<![=/<>])=(?![=>])|\.(?:eq|ne|lt|le|gt|ge|and|or|eqv|neqv|not)\.)\s*"
    )
    for m in re.finditer(pattern, code, re.IGNORECASE):
        if m.group(1):
            depth += 1
            edits.append((m.start(1) + 1, m.end(), ""))
            continue
        if m.group(2):
            depth -= 1
            edits.append((m.start(), m.start(2), ""))
            continue
        token = m.group(3).lower()
        # a blank before the first thing of the statement, or after its last, is not
        # part of the spacing
        start = m.start() if m.start() > 0 else m.start(3)
        end = m.end() if m.end() < len(code) else m.end(3)
        # an operator named in operator(...) has no blanks around it
        if re.search(r"\boperator\s*\(\s*$", code[: m.start(3)], re.IGNORECASE):
            edits.append((m.start(), m.end(), token))
            continue
        # inside brackets, :: follows the type of an allocate statement or an array
        # constructor, whereas an array section such as a(::2) is left as written
        type_spec = re.search(
            r"(\)|\b(integer|real|logical|complex|character|precision))\s*$",
            code[: m.start(3)],
            re.IGNORECASE,
        )
        if token == "::":
            if depth == 0 or type_spec:
                edits.append((start, end, " :: "))
        elif token == ",":
            edits.append((start, end, ", "))
        elif token in ("//", "%"):
            edits.append((start, end, token))
        elif token == ":":
            # tight in an array section or a bound, followed by one blank after only or
            # the name of a construct
            edits.append((start, end, ":" if depth > 0 else ": "))
        elif token == ".not.":
            # what precedes a unary operator is spaced by the token before it
            edits.append((m.start(3), end, ".not. "))
        elif token == "=":
            edits.append((start, end, "=" if depth > 0 else " = "))
        else:
            edits.append((start, end, f" {token} "))
    return edits


def format_literals(code: str) -> Set[int]:
    """
    this function returns where the string literals giving the format of an input/output
    statement start: the second positional or the fmt argument of read and write, and
    the first item of print and of read without brackets
    """
    starts = set()
    for start in statement_starts(code):
        m = re.match(r"(write|read)\s*\(", code[start:], re.IGNORECASE)
        if m:
            args = arguments(code, start + m.end() - 1)
            for k, (first, keyword) in enumerate(args):
                if keyword == "fmt":
                    starts.add(next_content(code, code.index("=", first) + 1))
                elif keyword is None and k == 1:
                    starts.add(first)
        m = re.match(r"(print|read)\s+(?=['\"])", code[start:], re.IGNORECASE)
        if m:
            starts.add(start + m.end())
    return starts


def requote(joined: Joined) -> Joined:
    """
    this function returns the statement with the format of an input/output statement
    written in single quotes and every other string literal in double quotes
    """
    formats = format_literals(joined.code)
    edits = []
    for first, last in literals(joined.text):
        quote = "'" if first in formats else '"'
        have = joined.text[first]
        if have != quote:
            inner = joined.text[first + 1 : last].replace(have * 2, have)
            edits.append(
                (first, last + 1, quote + inner.replace(quote, quote * 2) + quote)
            )
    if not edits:
        return joined
    edited = apply_edits(joined, edits)
    code = code_part(edited.text)
    if len(code) != len(edited.text):
        return joined
    return Joined(edited.text, code, edited.breaks, joined.reflowable)


def collapse_blanks(code: str) -> List[Tuple[int, int, str]]:
    """
    this function returns the replacements writing every run of blanks outside string
    literals as a single blank, so that code is never aligned with padding
    """
    spans = literals(code)
    edits = []
    for m in re.finditer(r"(?<=\S) {2,}(?=\S)", code):
        if not any(first < m.start() < last for first, last in spans):
            edits.append((m.start(), m.end(), " "))
    return edits


def normalize(joined: Joined) -> Joined:
    """
    this function returns the statement with its operators, punctuation and keywords
    written as rules 7 and 9 say
    """
    joined = requote(joined)
    passes = [
        space_operators,
        lower_keywords,
        keyword_words,
        keyword_brackets,
        space_punctuation,
        keyword_arguments,
        real_literals,
        unary_signs,
        collapse_blanks,
    ]
    for rule in passes:
        joined = apply_edits(joined, rule(joined.code))
    return joined


# the structure of a joined statement


DOTTED = re.compile(r"\.(and|or|eqv|neqv|eq|ne|lt|le|gt|ge)\.", re.IGNORECASE)

# how tightly a comma or binary operator binds its operands
PRECEDENCE = {",": 0, ".eqv.": 1, ".neqv.": 1, ".or.": 2, ".and.": 3, "//": 6}
PRECEDENCE.update({op: 5 for op in ("==", "/=", "<", "<=", ">", ">=")})
PRECEDENCE.update({f".{op}.": 5 for op in ("eq", "ne", "lt", "le", "gt", "ge")})
PRECEDENCE.update({"+": 7, "-": 7, "*": 8, "/": 8, "**": 9})

# the operators whose two operands are kept on one line under rule 4
COMPARISONS = {"==", "/=", "<", "<=", ">", ">=", ".eqv.", ".neqv."}
COMPARISONS.update(f".{op}." for op in ("eq", "ne", "lt", "le", "gt", "ge"))


@dataclass
class Structure:
    """
    this class holds the anchors, break candidates and bracket groups of a joined
    statement, which its layouts are chosen from
    """

    # position -> events at that position: ("open", content, start, close, lead, leaf,
    # array, single) for an opening bracket, ("close",) for a closing one, ("anchor",
    # content, kind) for any other anchor and ("entity_end",) for a comma ending an
    # entity of a declaration
    events: Dict[int, List[tuple]]
    candidates: List[int]  # positions a line may be broken at, ascending
    inner: Dict[int, int]  # breaks -> number of groups with a looser break in them
    comparisons: List[Tuple[int, int]]  # (start, end) of every comparison
    strings: Dict[int, str]  # breaks inside a string literal -> its quote character
    lengths: Dict[int, int]  # breaks inside a string literal -> length of the literal
    multi_entity: bool  # whether this is a declaration of several entities
    matrices: List[Tuple[int, int, List[int]]]  # (open, close, row ends)
    # (start, end) of every entity with an initializer in a declaration of several
    # entities, where end is the comma after it or the end of the statement
    entity_spans: List[Tuple[int, int]]


@dataclass
class Context:
    """
    this class holds what the declarations of a file tell about the layout of other
    statements: the values of integer named constants and the extents of arrays, each
    only for names whose declarations all agree
    """

    constants: Dict[str, int] = field(default_factory=dict)
    extents: Dict[str, List[int]] = field(default_factory=dict)


def integer_value(expr: str, context: Context) -> Optional[int]:
    """
    this function returns the value of an integer literal or named constant, or None
    when it cannot be resolved
    """
    expr = expr.strip().lower()
    m = re.fullmatch(r"(\d+)(_\w+)?", expr)
    return int(m.group(1)) if m else context.constants.get(expr)


def resolve_extents(code: str, at: int, context: Context) -> Optional[List[int]]:
    """
    this function returns the integers of an array constructor such as [n, 3] or the
    extents of shape(a) starting at a position, or None when they cannot be resolved
    """
    m = re.match(r"shape\s*\(\s*(\w+)\s*\)", code[at:], re.IGNORECASE)
    if m:
        return context.extents.get(m.group(1).lower())
    if code[at : at + 1] != "[":
        return None
    values = [integer_value(code[a:b], context) for a, b in pieces(code, at)]
    return None if None in values else values


def pieces(code: str, o: int) -> List[Tuple[int, int]]:
    """
    this function returns the start and end of every comma-separated piece of the
    bracket group opened at o
    """
    result, depth, start = [], 0, o + 1
    close = closing(code, o)
    for i in range(o + 1, close + 1):
        if i < close and code[i] in "([":
            depth += 1
        elif i < close and code[i] in ")]":
            depth -= 1
        elif i == close or (depth == 0 and code[i] == ","):
            result.append((start, i))
            start = i + 1
    return result


def file_context(codes: List[str]) -> Context:
    """
    this function returns the values of the integer named constants and the extents of
    the arrays declared in a file, leaving out names declared differently in two places
    """
    context, conflicts = Context(), set()

    def record(table: dict, name: str, value) -> None:
        if name in conflicts:
            return
        if value is None or table.get(name, value) != value:
            table.pop(name, None)
            conflicts.add(name)
        else:
            table[name] = value

    for code in codes:
        depth, decl = 0, None
        for i, ch in enumerate(code):
            if ch in "([":
                depth += 1
            elif ch in ")]":
                depth -= 1
            elif depth == 0 and code.startswith("::", i):
                decl = i
                break
        if decl is None or re.match(r"\s*use\b", code, re.IGNORECASE):
            continue
        attributes = code[:decl].lower()
        constant = re.match(r"\s*integer\b", attributes) and "parameter" in attributes
        listed = "(" + code[decl + 2 :] + ")"
        for a, b in pieces(listed, 0):
            entity = listed[a:b].strip()
            m = re.match(r"(\w+)\s*(\(([^=]*)\))?\s*(=(?!>)\s*(.*))?$", entity)
            if not m:
                continue
            name = m.group(1).lower()
            if m.group(2):
                bounds = "(" + m.group(3) + ")"
                extents = [
                    integer_value(bounds[c:d], context) for c, d in pieces(bounds, 0)
                ]
                record(context.extents, name, None if None in extents else extents)
            if constant and m.group(4) and not m.group(2):
                record(context.constants, name, integer_value(m.group(5), context))
    return context


def structure(code: str, text: str, context: Context) -> Structure:
    """
    this function finds the anchors and break candidates of a joined statement
    """
    events: Dict[int, List[tuple]] = {}
    candidates: Set[int] = set()

    def add(i: int, event: tuple) -> None:
        events.setdefault(i, []).append(event)

    lower = code.lower()
    end = len(code.rstrip())
    is_use = re.match(r"\s*use\b", lower) is not None

    # matching brackets
    match, stack = {}, []
    for i, ch in enumerate(code):
        if ch in "([":
            stack.append(i)
        elif ch in ")]" and stack:
            match[stack.pop()] = i

    # a declaration has a :: outside brackets
    depth, decl_at = 0, None
    for i, ch in enumerate(code):
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        elif depth == 0 and code.startswith("::", i) and not is_use:
            decl_at = i
            break

    # parentheses that only group an expression, as opposed to those of a call, an
    # array reference or a statement keyword such as if
    def grouping(o: int) -> bool:
        j = o - 1
        while j >= 0 and code[j] == " ":
            j -= 1
        return code[o] == "(" and (j < 0 or not (code[j].isalnum() or code[j] in "_%"))

    # how tightly the operator or comma every break follows binds, in the bracket group
    # it is in and in the groups around it, where it binds tighter than any operator or
    # comma outside its bracket group
    binds: Dict[int, List[Tuple[int, int]]] = {}
    opening: Dict[int, List[Tuple[int, int]]] = {}
    groups = [-1]

    def bind(rank: int) -> List[Tuple[int, int]]:
        entries = [(rank, groups[-1])]
        for k in range(len(groups) - 1, 0, -1):
            rank += 10
            entries.append((rank, groups[k - 1]))
        return entries

    depth, assigned, in_list, entities = 0, None, False, 0
    initializers, entity_starts, entity_ends = [], [], []
    # every operator, comma and assignment as (start, end, rank, group, operator)
    ops: List[Tuple[int, int, int, int, str]] = []
    i = 0
    while i < end:
        ch = code[i]
        if ch in "([":
            groups.append(i)
            # the group a bracket belongs to starts at the name before it
            start = i
            while start > 0 and (code[start - 1].isalnum() or code[start - 1] in "_%"):
                start -= 1
            add(i, ("open", next_content(code, i + 1), start, match.get(i, end)))
            depth += 1
            candidates.add(i + 1)
            # a break directly after the bracket is inside its group, and so binds
            # tighter than anything in the groups around it
            opening[i + 1] = bind(PRECEDENCE[","])[1:]
            i += 1
            continue
        if ch in ")]":
            add(i, ("close",))
            depth -= 1
            if len(groups) > 1:
                groups.pop()
            i += 1
            continue
        if ch == ",":
            candidates.add(i + 1)
            binds[i + 1] = bind(PRECEDENCE[","])
            ops.append((i, i + 1, PRECEDENCE[","], groups[-1], ","))
            if in_list and depth == 0:
                add(i, ("entity_end",))
                entity_ends.append(i)
                entity_starts.append(next_content(code, i + 1))
                entities += 1
            i += 1
            continue
        if decl_at is not None and i == decl_at:
            add(i + 1, ("anchor", next_content(code, i + 2), "list"))
            candidates.add(i + 2)
            in_list, entities = True, 1
            entity_starts.append(next_content(code, i + 2))
            i += 2
            continue
        if is_use and depth == 0 and lower.startswith("only", i):
            colon = next_content(code, i + 4)
            if code[colon : colon + 1] == ":":
                add(colon, ("anchor", next_content(code, colon + 1), "list"))
                candidates.add(colon + 1)
                i = colon + 1
                continue
        if ch == "=" and depth == 0 and not is_use:
            prev = code[i - 1] if i > 0 else " "
            nxt = code[i + 1] if i + 1 < len(code) else " "
            if nxt != "=" and prev not in "=/<>":
                op_end = i + 2 if nxt == ">" else i + 1
                ops.append((i, op_end, -1, groups[-1], "="))
                if in_list:
                    add(op_end - 1, ("anchor", next_content(code, op_end), "init"))
                    candidates.add(op_end)
                    initializers.append(op_end)
                elif not assigned and decl_at is None:
                    add(op_end - 1, ("anchor", next_content(code, op_end), "assign"))
                    candidates.add(op_end)
                    assigned = op_end
                i = op_end
                continue
        # binary operators
        op = None
        m = DOTTED.match(code, i)
        if m:
            op = m.end()
        elif code.startswith(("//", "**", "==", "/=", "<=", ">="), i):
            if operand_before(code, i):
                op = i + 2
        elif ch in "+-*/<>" and operand_before(code, i):
            if not (ch == "/" and code[i + 1 : i + 2] == ")"):
                op = i + 1
        if op is not None:
            candidates.add(op)
            name = lower[i:op].strip()
            binds[op] = bind(PRECEDENCE[name])
            ops.append((i, op, PRECEDENCE[name], groups[-1], name))
            i = op
            continue
        i += 1

    # in a declaration of several entities, the initializer of an entity binds tighter
    # than the commas separating the entities
    if entities > 1:
        for b in initializers:
            binds[b] = [(1, -1)]

    # a break after the condition of a single-line if or where statement, which is
    # looser than any break in the statement it controls
    m = re.match(r"\s*(if|where)\s*\(", lower)
    if m and match.get(m.end() - 1) is not None:
        after = match[m.end() - 1] + 1
        if not re.match(r"\s*then\s*$", lower[after:]):
            candidates.add(after)
            binds[after] = [(-2, -1)]
            # the assignment it controls, which binds tighter than the condition
            if assigned is not None and assigned > after:
                binds[assigned] = [(-1, -1)]

    # breaks inside string literals, after a blank between two words, where the word
    # after the break is not made of punctuation only, such as | or =
    strings: Dict[int, str] = {}
    lengths: Dict[int, int] = {}
    for first, last in literals(text):
        if aligned_by_hand(text, first, last):
            continue
        for b in range(first + 3, last):
            if text[b - 1] != " " or text[b] == " " or text[b - 2] == " ":
                continue
            word = re.match(r"[^ ]*", text[b:last]).group(0)
            if re.search(r"[^\W_]", word):
                strings[b] = text[first]
                lengths[b] = last - first + 1

    # a break before a lone then, which is looser than any break in the condition, and
    # before the result or bind of a procedure
    m = re.search(r"\)\s*(then)\s*$", lower)
    if m:
        candidates.add(m.start(1))
        binds[m.start(1)] = [(-2, -1)]
    for kw in ("result", "bind"):
        for m in re.finditer(rf"\)\s*({kw})\s*\(", lower):
            candidates.add(m.start(1))

    # a break that splits an operand of a looser operator or comma in the same group
    loosest: Dict[int, int] = {}
    for entries in binds.values():
        for rank, group in entries:
            loosest[group] = min(rank, loosest.get(group, rank))
    inner = {
        b: sum(rank > loosest.get(group, rank) for rank, group in entries)
        for b, entries in list(binds.items()) + list(opening.items())
    }

    candidates = {b for b in candidates if 0 < b < end and next_content(code, b) < end}
    candidates.update(strings)

    # the extent of every comparison, bounded on either side by the nearest operator,
    # comma or assignment of its bracket group that binds more loosely, or else by the
    # group itself
    extents = []
    for first, last, rank, group, name in ops:
        if name not in COMPARISONS:
            continue
        left = group + 1 if group >= 0 else 0
        right = match.get(group, end) if group >= 0 else end
        for other_first, other_last, other_rank, other_group, _ in ops:
            if other_group != group or other_rank >= rank:
                continue
            if other_last <= first:
                left = max(left, other_last)
            elif other_first >= last:
                right = min(right, other_first)
        extents.append((next_content(code, left), len(code[:right].rstrip())))

    # an array constructor, an array section, or the shape of an array in a declaration
    # or an allocate statement
    allocating = re.match(r"\s*(de)?allocate\s*\(", lower)

    def array_like(o: int, c: int, enclosing: List[int]) -> bool:
        if code[o] == "[":
            return True
        depth = 0
        for j in range(o + 1, c):
            if code[j] in "([":
                depth += 1
            elif code[j] in ")]":
                depth -= 1
            elif depth == 0 and code[j] == ":" and ":" not in code[j - 1 : j + 2 : 2]:
                return True
        if decl_at is not None and o > decl_at and not enclosing:
            return True
        return bool(allocating) and len(enclosing) == 1

    # a bracket group holding a single argument, index or expression
    def single(o: int, c: int) -> bool:
        depth = 0
        for j in range(o + 1, c):
            if code[j] in "([":
                depth += 1
            elif code[j] in ")]":
                depth -= 1
            elif depth == 0 and code[j] == ",":
                return False
        return True

    # a bracket group can be moved to the next line whole, together with what joins it
    # to the last possible break before it (such as .not.), when that break lies inside
    # the bracket group enclosing it
    ordered = sorted(candidates - set(strings))
    for i, evs in events.items():
        for k, ev in enumerate(evs):
            if ev[0] == "open":
                before = [b for b in ordered if next_content(code, b) <= ev[2]]
                lead = next_content(code, before[-1]) if before else None
                enclosing = [o for o, c in match.items() if o < ev[2] and c > ev[3]]
                if lead is not None and enclosing and lead <= max(enclosing):
                    lead = None
                leaf = not any(i < o < ev[3] and not grouping(o) for o in match)
                array = array_like(i, ev[3], enclosing)
                evs[k] = ev + (lead, leaf, array, single(i, ev[3]))

    # array constructors passed to reshape, laid out one matrix row per line, where a
    # row holds as many items as the extent filled first
    matrices = []
    for m in re.finditer(r"\breshape\s*\(", lower):
        args = arguments(code, m.end() - 1)
        open_at = args[0][0] if args else None
        if open_at is None or code[open_at] != "[" or open_at not in match:
            continue
        close_at = match[open_at]
        depth, commas = 0, []
        for i in range(open_at + 1, close_at):
            if code[i] in "([":
                depth += 1
            elif code[i] in ")]":
                depth -= 1
            elif code[i] == "," and depth == 0:
                commas.append(i)
        n_items = len(commas) + 1
        shape = order = None
        for k, (first, keyword) in enumerate(args):
            value = first
            if keyword is not None:
                value = next_content(code, code.index("=", first) + 1)
            if keyword == "shape" or (keyword is None and k == 1):
                shape = resolve_extents(code, value, context)
            elif keyword == "order" or (keyword is None and k == 3):
                order = resolve_extents(code, value, context)
        row = None
        if shape and not order:
            row = shape[0]
        elif shape and order and 1 <= order[0] <= len(shape):
            row = shape[order[0] - 1]
        if row and 1 < row < n_items and n_items % row == 0:
            row_ends = [commas[k * row - 1] + 1 for k in range(1, n_items // row)]
            matrices.append((open_at, close_at, row_ends))

    return Structure(
        events,
        sorted(candidates),
        inner,
        extents,
        strings,
        lengths,
        entities > 1,
        matrices,
        [
            (start, stop)
            for start, stop in zip(entity_starts, entity_ends + [end])
            if entities > 1 and any(start < b <= stop for b in initializers)
        ],
    )


# the layout search


# flags of an open bracket group under rule 4: it fits whole on a line, contains no
# other bracket group apart from grouping parentheses, is array-like, or holds a single
# argument or index
FITS, LEAF, ARRAY, SINGLE = 1, 2, 4, 8

# an open anchor: its kind, the column of the first thing after it (None when that
# hangs), its hanging indentation, its rule 4 flags and the position of the first thing
# after it
Anchor = Tuple[str, Optional[int], int, int, int]

# the cost of a layout in the order it is minimized (array-like groups and comparisons
# split, lines, other groups split, operands split), and its breaks, negated so that
# the latest breaks sort first
Score = Tuple[int, int, int, int, Tuple[int, ...]]


class Layout:
    """
    this class holds the layouts of one joined statement starting at a given column
    """

    def __init__(
        self,
        joined: Joined,
        base: int,
        s: Structure,
        matrix: bool = False,
        entity: bool = False,
    ):
        self.j = joined
        self.base = base
        self.s = s
        self.end = len(joined.code.rstrip())
        self.forced: List[int] = []
        self.inside: List[Tuple[int, int]] = []
        self.rows = matrix and bool(self.s.matrices)
        if self.rows:
            for open_at, close_at, row_ends in self.s.matrices:
                self.forced += row_ends
                self.inside.append((open_at, close_at))
        # an entity with an initializer that fits whole on a continuation line is kept
        # whole like a comparison, whereas one that does not starts a new line, as does
        # the entity after it
        fitting, long = [], []
        for a, c in self.s.entity_spans:
            fits = base + HANG + c - a + len(", &") <= LIMIT
            (fitting if fits else long).append((a, c))
        self.entities = entity and bool(long)
        if self.entities:
            for a, c in long:
                before = len(joined.code[:a].rstrip())
                for b in (before, c + 1):
                    if b in self.s.candidates and b not in self.forced:
                        self.forced.append(b)
        self.forced.sort()
        # the start of the bracket group, name included, whose contents begin at a
        # position
        self.starts = {
            ev[1]: ev[2]
            for evs in self.s.events.values()
            for ev in evs
            if ev[0] == "open"
        }
        # a string literal that fits whole on a continuation line is split only when
        # that saves a line
        self.whole_strings = {
            b
            for b, length in self.s.lengths.items()
            if base + HANG + length + len(", &") <= LIMIT
        }
        self.candidates = [
            b
            for b in self.s.candidates
            if b in self.forced or not any(a < b < c for a, c in self.inside)
        ]
        # the comparisons a break splits, counting only those that fit whole on a
        # continuation line with a hanging indent
        whole_comparisons = [
            (a, c)
            for a, c in self.s.comparisons
            if base + HANG + c - a + len(" &") <= LIMIT
        ] + fitting
        self.compared = {
            b: sum(a < b < c for a, c in whole_comparisons) for b in self.candidates
        }
        # whether the line after a break starts with then, result or bind
        self.keyword_starts: Dict[int, bool] = {}

    def scan(
        self, s: int, b: int, indent: int, stack: Tuple[Anchor, ...]
    ) -> Tuple[Anchor, ...]:
        """
        this function returns the open anchors after laying out positions s to b on a
        line starting at the given indentation
        """

        # a line continuing a string literal starts with its quote
        pre = int(s in self.s.strings)

        def col(i: int) -> int:
            return indent + pre + i - s

        anchors = list(stack)
        for i in range(s, b):
            for ev in self.s.events.get(i, ()):
                if ev[0] == "open":
                    _, content, start, close, lead, leaf, array, one = ev
                    fits = False
                    # a group holding matrix rows is split by them
                    rows = any(i <= a and c <= close for a, c in self.inside)
                    if lead is not None and lead >= s and not rows:
                        # where the group would start if moved to the next line
                        if lead == s:
                            at = indent
                        elif anchors:
                            _, anchor, hang, _, top_content = anchors[-1]
                            if anchor is None or lead <= top_content:
                                at = hang
                            else:
                                at = anchor
                        else:
                            at = self.base + HANG
                        # whether it fits there or on a continuation line with a
                        # hanging indent, which an earlier break can give it
                        at = min(at, self.base + HANG)
                        tail = 2 if close + 1 < self.end else 0
                        fits = at + close - lead + 1 + tail <= LIMIT
                    # an innermost group that fits whole where it stands is not split
                    # either
                    tail = 2 if close + 1 < self.end else 0
                    here = col(start) + close - start + 1 + tail <= LIMIT
                    fits = fits or (leaf and not rows and here)
                    whole = 0
                    if fits:
                        whole = FITS | LEAF * leaf | ARRAY * array | SINGLE * one
                    anchor = col(content) if content < b else None
                    if anchor is None:
                        # a bracket whose first argument hangs hangs with it
                        pos = start
                        for k in range(len(anchors) - 1, -1, -1):
                            kind, _, hang, outer_whole, outer_content = anchors[k]
                            if kind != "bracket" or outer_content != pos:
                                break
                            anchors[k] = (kind, None, hang, outer_whole, outer_content)
                            pos = self.starts[outer_content]
                    anchors.append(("bracket", anchor, indent + HANG, whole, content))
                elif ev[0] == "close":
                    if anchors and anchors[-1][0] == "bracket":
                        anchors.pop()
                elif ev[0] == "anchor":
                    _, content, kind = ev
                    names = next((a for a in anchors if a[0] == "list"), None)
                    if kind == "init" and self.s.multi_entity and names is not None:
                        hang = (names[2] if names[1] is None else names[1]) + HANG
                    else:
                        hang = self.base + HANG
                    anchor = col(content) if content < b else None
                    anchors.append((kind, anchor, hang, 0, content))
                elif ev[0] == "entity_end":
                    if anchors and anchors[-1][0] == "init":
                        anchors.pop()
        return tuple(anchors)

    def next_indent(self, b: int, stack: Tuple[Anchor, ...]) -> int:
        """
        this function returns the indentation of the line following a break
        """
        # a line starting with then, result or bind hangs from the statement
        if b not in self.keyword_starts:
            rest = self.j.code[next_content(self.j.code, b) :]
            keyword = re.match(r"then\s*$|(result|bind)\s*\(", rest, re.IGNORECASE)
            self.keyword_starts[b] = keyword is not None
        if self.keyword_starts[b]:
            return self.base + HANG
        if stack:
            _, anchor, hang, _, _ = stack[-1]
            return hang if anchor is None else anchor
        return self.base + HANG

    def best(self) -> Optional[List[int]]:
        """
        this function returns the break positions of the best layout, or None when no
        layout keeps within the column limit
        """
        text, end = self.j.text, self.end
        first = next_content(self.j.code, 0)

        # a single line is the best layout wherever it fits, since every other one has
        # more lines and none splits fewer groups
        if not self.forced and self.base + len(text[first:end].rstrip()) <= LIMIT:
            return []

        @lru_cache(maxsize=None)
        def search(s: int, indent: int, stack: Tuple[Anchor, ...]) -> Optional[Score]:
            options: List[Score] = []
            forced = next((f for f in self.forced if f > s), None)
            start = indent + int(s in self.s.strings)
            if forced is None and start + len(text[s:end].rstrip()) <= LIMIT:
                options.append((0, 1, 0, 0, ()))
            for b in reversed(self.candidates):
                if b <= s or (forced is not None and b > forced):
                    continue
                if b in self.s.strings:
                    # the literal is closed and concatenated with its continuation
                    width = start + len(text[s:b]) + len('"// &')
                else:
                    width = start + len(text[s:b].rstrip()) + len(" &")
                if width > LIMIT:
                    continue
                new_stack = self.scan(s, b, indent, stack)
                # splitting an array-like group that fits whole, and splitting an
                # innermost group or hanging the contents of any group that fits whole
                arrays = groups = 0
                for _, anchor, _, whole, _ in new_stack:
                    hung = anchor is None
                    if whole & ARRAY or (whole & SINGLE and hung):
                        arrays += 1
                    elif whole & LEAF or (whole and hung):
                        groups += 1
                split = self.s.inner.get(b, 0)
                arrays += self.compared.get(b, 0)
                groups += b in self.whole_strings
                n = b if b in self.s.strings else next_content(self.j.code, b)
                sub = search(n, self.next_indent(b, new_stack), new_stack)
                if sub is not None:
                    options.append(
                        (
                            sub[0] + arrays,
                            sub[1] + 1,
                            sub[2] + groups,
                            sub[3] + split,
                            (-b,) + sub[4],
                        )
                    )
            # fewest array-like groups and comparisons split, then fewest lines, then
            # fewest other groups split against rule 4, then fewest split operands, then
            # the latest first break, the latest second break, ...
            return min(options) if options else None

        result = search(first, self.base, ())
        return None if result is None else [-b for b in result[4]]

    def render(self, breaks: List[int]) -> List[str]:
        """
        this function returns the lines of the statement broken at the given positions
        """
        lines, indent = [], self.base
        stack: Tuple[Anchor, ...] = ()
        s = next_content(self.j.code, 0)
        for b in breaks:
            line = " " * indent + self.s.strings.get(s, "")
            if b in self.s.strings:
                line += self.j.text[s:b] + self.s.strings[b] + "// &"
            else:
                line += self.j.text[s:b].rstrip() + " &"
            lines.append(line)
            stack = self.scan(s, b, indent, stack)
            indent = self.next_indent(b, stack)
            s = b if b in self.s.strings else next_content(self.j.code, b)
        start = " " * indent + self.s.strings.get(s, "")
        lines.append(start + self.j.text[s : self.end].rstrip())
        return lines


def best_layout(joined: Joined, base: int, context: Context) -> Optional[List[str]]:
    """
    this function returns the best layout of a statement, dropping the matrix rows and
    then the entities on their own lines when no layout keeps them, and keeping string
    literals as they are written when merging them leaves no layout at all
    """
    for statement in (merge_literals(requote(joined)), joined):
        statement = normalize(statement)
        s = structure(statement.code, statement.text, context)
        for matrix, entity in (
            (True, True),
            (True, False),
            (False, True),
            (False, False),
        ):
            layout = Layout(statement, base, s, matrix, entity)
            if (matrix and not layout.rows) or (entity and not layout.entities):
                continue
            breaks = layout.best()
            if breaks is not None:
                return layout.render(breaks)
    return None


# checking and fixing


@dataclass
class Issue:
    """
    this class holds an issue of a statement, together with what fixes it
    """

    line: int
    kind: str
    message: str
    stmt: Statement
    replacement: Optional[List[str]] = None  # new lines of the whole statement
    indent: Optional[int] = None  # the column the line should start at


def in_scope(stmt: Statement, touched: Optional[Set[int]]) -> bool:
    """
    this function returns whether a statement is to be checked, which it is when it
    contains a touched line or when every line counts as touched
    """
    return touched is None or bool(
        touched.intersection(range(stmt.first, stmt.last + 1))
    )


def analyse(
    lines: List[str], touched: Optional[Set[int]], layout: bool = True
) -> List[Issue]:
    """
    this function returns the issues of a file, restricted to statements containing a
    touched line when a set of touched lines is given, and leaving out the issues within
    statements, whose layout search is what takes time, unless asked for
    """
    issues = []
    stmts = statements(lines)
    joins = [join_statement(lines, stmt) for stmt in stmts]
    codes = [joined.code for joined in joins]
    bases = expected_indents(codes)
    ends, one, none = block_rules(lines, stmts, codes)
    context = file_context(codes)

    # no line ends in blanks, no file starts or ends with a blank line, no two blank
    # lines follow each other, and the file ends with a newline; what follows the final
    # newline is no line
    terminated = bool(lines) and not lines[-1].strip()
    last = len(lines) - 1 if terminated else len(lines)
    if lines and not terminated:
        here = Statement(last - 1, last - 1, [])
        issues.append(Issue(last - 1, "newline", "no newline at end of file", here))
    content = [i for i in range(last) if lines[i].strip()]
    superfluous = set()
    for i in range(last):
        if not lines[i].strip():
            outside = not content or i < content[0] or i > content[-1]
            if outside or not lines[i - 1].strip():
                superfluous.add(i)
    for i in none:
        j = i + 1
        while j < last and not lines[j].strip():
            superfluous.add(j)
            j += 1
    for i in range(len(lines)):
        line, here = lines[i], Statement(i, i, [])
        if in_scope(here, touched) and line != line.rstrip():
            issues.append(Issue(i, "whitespace", "trailing blanks", here))
        if in_scope(here, touched) and i in superfluous:
            issues.append(Issue(i, "blank", "superfluous blank line", here))
        gap = Statement(i, i + 1, [])
        if (
            i in one
            and i + 1 < last
            and lines[i + 1].strip()
            and in_scope(gap, touched)
        ):
            issues.append(Issue(i, "gap", "missing blank line after this one", here))

    # comment lines between statements start where the statement after them does, while
    # those inside a continued statement move with it
    inside = {i for stmt in stmts for i in range(stmt.first, stmt.last + 1)}
    k = 0
    for i, line in enumerate(lines):
        if not line.lstrip().startswith("!"):
            continue
        while k < len(stmts) and stmts[k].first < i:
            k += 1
        column = bases[k] if k < len(stmts) else 0
        if i in inside:
            column = indent_of(line)
        comment = Statement(i, i, [])
        want = " " * column + comment_text(line.strip())
        if want != line.rstrip() and in_scope(comment, touched):
            message = f"comment starts at {indent_of(line)}, expected {column}"
            if indent_of(line) == column:
                message = "comment text not spaced as rule 10 says"
            issues.append(Issue(i, "comment", message, comment, [want]))

    for n, (stmt, base) in enumerate(zip(stmts, bases)):
        if not in_scope(stmt, touched):
            continue
        have = indent_of(lines[stmt.first])
        if have != base:
            message = f"starts at {have}, expected {base}"
            issues.append(Issue(stmt.first, "indent", message, stmt, indent=base))
            continue
        if not layout:
            continue
        joined = joins[n]
        for literal in kindless_reals(joined.code):
            message = f"real literal {literal} needs a kind, fix by hand"
            issues.append(Issue(stmt.first, "literal", message, stmt))
        renamed = n in ends and joined.code.strip() != ends[n]
        if renamed:
            start = next_content(joined.code, 0)
            joined = apply_edits(joined, [(start, len(joined.code), ends[n])])
        if not joined.reflowable:
            # keep the breaks and the comments, only fix where the lines start and how
            # they are spaced
            normal = normalize(joined)
            s = structure(normal.code, normal.text, context)
            want = Layout(normal, base, s).render(normal.breaks)
            current = [line.rstrip() for line in lines[stmt.first : stmt.last + 1]]
            replacement = list(current)
            for i, want_line in zip(stmt.code, want):
                comment = lines[i][len(code_part(lines[i])) :].strip()
                if comment:
                    comment = "  " + comment_text(comment)
                replacement[i - stmt.first] = want_line + comment
            for i in stmt.code:
                line = replacement[i - stmt.first]
                if len(line) > LIMIT:
                    message = f"{len(line)} columns, statement not reflowable"
                    issues.append(Issue(i, "overlong", message, stmt))
            if replacement != current:
                same_text = [line.strip() for line in replacement] == [
                    line.strip() for line in current
                ]
                kind = "name" if renamed else "align" if same_text else "spacing"
                message = f"{kind}, statement not reflowable"
                issues.append(Issue(stmt.first, kind, message, stmt, replacement))
            continue
        current = [lines[i].rstrip() for i in range(stmt.first, stmt.last + 1)]
        best = best_layout(joined, base, context)
        if best is None:
            for i in stmt.code:
                if len(lines[i].rstrip()) > LIMIT:
                    message = f"{len(lines[i].rstrip())} columns, no layout fits"
                    issues.append(Issue(i, "overlong", message, stmt))
            continue
        if best != current:
            if renamed:
                kind = "name"
            elif any(len(line) > LIMIT for line in current):
                kind = "overlong"
            elif [line.strip() for line in best] == [line.strip() for line in current]:
                kind = "align"
            elif [re.sub(r"\s", "", line).lower() for line in best] == [
                re.sub(r"\s", "", line).lower() for line in current
            ]:
                kind = "spacing"
            else:
                kind = "layout"
            message = f"{len(current)} line(s), best layout has {len(best)}"
            issues.append(Issue(stmt.first, kind, message, stmt, best))
    return issues


def fix(
    lines: List[str], touched: Optional[Set[int]]
) -> Tuple[List[str], int, Optional[List[Issue]]]:
    """
    this function returns the lines with every fixable issue fixed, the number of fixes
    and, when there was nothing to fix, the issues that remain; blanks and blank lines
    are fixed first, since removing lines moves the statements, and then comments and
    the indentation of statements, alone, since the layout of a statement depends on the
    column it starts at
    """
    issues = analyse(lines, touched, layout=False)
    blank = {issue.line for issue in issues if issue.kind == "blank"}
    stripped = {issue.line for issue in issues if issue.kind == "whitespace"}
    gaps = {issue.line for issue in issues if issue.kind == "gap"}
    newline = any(issue.kind == "newline" for issue in issues)
    if blank or stripped or gaps or newline:
        result = []
        for i, line in enumerate(lines):
            if i not in blank:
                result.append(line.rstrip() if i in stripped else line)
            if i in gaps:
                result.append("")
        return result + [""] * newline, len(blank | stripped | gaps) + newline, None
    changed = 0
    for issue in issues:
        if issue.kind == "comment":
            lines[issue.line] = issue.replacement[0]
            changed += 1
            continue
        if issue.kind != "indent":
            continue
        stmt = issue.stmt
        shift = issue.indent - indent_of(lines[stmt.first])
        for i in range(stmt.first, stmt.last + 1):
            if lines[i].strip():
                line = lines[i]
                lines[i] = " " * max(indent_of(line) + shift, 0) + line.lstrip()
        changed += 1
    if changed:
        return lines, changed, None
    issues = analyse(lines, touched)
    for issue in reversed(issues):
        stmt = issue.stmt
        if issue.replacement is not None:
            lines = lines[: stmt.first] + issue.replacement + lines[stmt.last + 1 :]
            changed += 1
    return lines, changed, None if changed else issues


# command line


def git(*args: str) -> str:
    """
    this function returns the output of a git command, stopping with its error when it
    fails, so that a failing command is never mistaken for one reporting no changes
    """
    result = subprocess.run(
        ["git", "--no-pager", *args], capture_output=True, text=True
    )
    if result.returncode != 0:
        # the first line holds the reason, what follows is usage help
        reason = (result.stderr.strip().splitlines() or ["unknown error"])[0]
        sys.exit(f"git {' '.join(args)} failed: {reason}")
    return result.stdout


def touched_lines(path: Path, ref: str) -> Optional[Set[int]]:
    """
    this function returns the 0-based indices of lines added or changed relative to a
    git reference, or None when the file is not tracked, in which case every line
    counts as touched
    """
    tracked = subprocess.run(
        ["git", "ls-files", "--error-unmatch", str(path)], capture_output=True
    )
    if tracked.returncode != 0:
        return None
    diff = git("diff", "--no-color", "--no-ext-diff", "-U0", ref, "--", str(path))
    result: Set[int] = set()
    for hunk in re.finditer(r"^@@ -\S+ \+(\d+)(?:,(\d+))? @@", diff, re.MULTILINE):
        start, count = int(hunk.group(1)), int(hunk.group(2) or 1)
        if count == 0:
            # a hunk that only deletes lines touches the lines on either side of it
            result.update(i for i in (start - 1, start) if i >= 0)
        else:
            result.update(range(start - 1, start - 1 + count))
    return result


def main() -> int:
    """
    this function checks the files named or found through git, fixes what is reported
    when asked to, and returns 1 when issues remain and 0 otherwise
    """
    parser = argparse.ArgumentParser(description=__doc__.strip().split("\n")[0])
    parser.add_argument("files", nargs="*", type=Path, help="Fortran files to check")
    parser.add_argument("--ref", default="HEAD", help="git reference (default HEAD)")
    parser.add_argument(
        "--all", action="store_true", help="check whole files, not only changes"
    )
    parser.add_argument("--fix", action="store_true", help="fix what is reported")
    parser.add_argument(
        "--show",
        action="store_true",
        help="print the best layout of every statement reported",
    )
    args = parser.parse_args()

    if args.files:
        files = args.files
    elif args.all and args.fix:
        parser.error("--all --fix needs the files to fix to be named explicitly")
    else:
        listings = [["ls-files", "-z"]]
        if not args.all:
            # changed files, and files git does not track yet
            listings = [
                ["diff", "--no-color", "--name-only", "--relative", "-z", args.ref],
                ["ls-files", "-z", "--others", "--exclude-standard"],
            ]
        out = []
        for listing in listings:
            out += git(*listing).split("\0")
        files = [Path(f) for f in out if f.endswith(".f90") and Path(f).exists()]

    remaining = 0
    for path in files:
        lines = path.read_text().split("\n")

        def scope() -> Optional[Set[int]]:
            return None if args.all else touched_lines(path, args.ref)

        issues = None
        if args.fix:
            total, written = 0, "\n".join(lines)
            for _ in range(10):
                lines, changed, issues = fix(lines, scope())
                if not changed:
                    break
                # never overwrite a change saved to the file in the meantime
                if path.read_text() != written:
                    print(f"{path}: changed while fixing, left as it is")
                    lines = path.read_text().split("\n")
                    break
                written = "\n".join(lines)
                path.write_text(written)
                total += changed
            else:
                print(f"{path}: fixes did not settle")
            if total:
                print(f"{path}: {total} fix(es) applied")
        if issues is None:
            issues = analyse(lines, scope())
        for issue in issues:
            print(f"{path}:{issue.line + 1}: {issue.kind}: {issue.message}")
            if args.show and issue.replacement:
                print("\n".join("    | " + line for line in issue.replacement))
            remaining += 1
    return 1 if remaining else 0


if __name__ == "__main__":
    sys.exit(main())
