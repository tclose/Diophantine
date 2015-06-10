
;;; Copyright (C) 2006 Aarre Laakso <aarre@pair.com>
;;;  
;;; This program is free software; you can redistribute it and/or modify
;;; it under the terms of the GNU General Public License as published by
;;; the Free Software Foundation; either version 2 of the License, or
;;; (at your option) any later version.
;;;  
;;; This program is distributed in the hope that it will be useful,
;;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;;; GNU General Public License for more details.
;;;  
;;; You should have received a copy of the GNU General Public License
;;; along with this program; if not, write to the Free Software
;;; Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
;;;  

;;; $Id: php2py.el,v 1.2 2006/11/23 12:09:06 aarre Exp $

;; Convert PHP code to Python
;;
;; NAME
;;     php2py.el
;;
;; USAGE
;;     M-x load-file php2py.el
;;     M-x php2py
;;
;; DESCRIPTION
;;     Performs the mechanical work of converting PHP code to Python.
;;
;;
;; CREDIT
;;
;;     Started from the following list of regular expressions
;;     published by Kieth Devens at
;;     http://keithdevens.com/weblog/archive/2003/Nov/22/convert-PHP-to-Python
;;
;;     remove "<?php" or "?>" anywhere they appear
;;     s/;;     $// #remove trailing semicolons
;;     s/\$// #remove dollar signs
;;     s/.=/+=/ #change concatenation to Python syntax
;;     s/array\(/[/ #though, that's only right in some cases
;;     s/list\(/\(/
;;     s/\s*=>/:/
;;     s/->/./
;;     s/if\((.*)\){/if $1:/ #remove any parentheses around ifs
;;     s/}// # Smiley (big smile)
;;     s/elseif/elif/
;;     s/else{/else:/
;;     s/echo/print/
;;     s/NULL/None/i
;;     s/++/+=1/
;;     s/is_numeric\(([^ ]*)\)/$1.isdigit()/
;;     s/strlen/len/
;;

;; rx.el is distributed with GNU Emacs. It seems to work just fine with XEmacs.
(require 'rx)

(defun php2py () (interactive)

  ;; remove "<?php"
  (goto-char (point-min))
  (query-replace-regexp (rx "<?php") "")

  ;; remove "?>"
  (goto-char (point-min))
  (query-replace-regexp (rx "?>") "")

  ;; remove trailing semicolons
  (goto-char (point-min))
  (perform-replace ";" "" t nil nil)

  ;; remove "$"
  (goto-char (point-min))
  (perform-replace "$" "" t nil nil)

  ;; change concatenation syntax
  (goto-char (point-min))
  (query-replace-regexp
   (rx "\" . ")
   "\" + ")
  (goto-char (point-min))
  (query-replace-regexp
   (rx " . \"")
   " + \"")
  (goto-char (point-min))
  (perform-replace ".=" "+=" t nil nil)

  ;; Change ++ to += 1
  (goto-char (point-min))
  (query-replace-regexp 
   (rx "++")
   " += 1")

  ;; Change && to and
  (goto-char (point-min))
  (query-replace-regexp
   (rx "&&")
   "and")

  ;; Change ! to not
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq (char "!")
	    (submatch (not-char "="))))
   " not \\1")

  ;; Add line continuations
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq
	(in ".+")
	(zero-or-more " ")
	line-end)) "+ \\")
  
  ;; change attribute access syntax
  (goto-char (point-min))
  (perform-replace "->" "." t nil nil)

  ;; Change block comments to docstrings
  (goto-char (point-min))
  (query-replace-regexp 
   (rx
    (seq "/*"
	 (submatch (minimal-match (zero-or-more anything)))
	 "*/"))
   "\"\"\" \\1 \"\"\"")

  ;; Remove extra asterisks from within block comments / docstrings
  (goto-char (point-min))
  (query-replace-regexp "*" "")

  ;; Change // comments to #
  (goto-char (point-min))
  (query-replace-regexp "//" "#")

  ;; change incude_once(...) to import ...
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "include_once("
	    (submatch (minimal-match (zero-or-more anything)))
	    ")"))
   "import \\1")

  ;; Change for (x=n; x<=m; x+=1) to for x in range(n,m+1)
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "for"
	    (zero-or-more " ")
	    (char "(")
	    (zero-or-more " ")
	    (submatch (one-or-more (not (in " ")))) ; x
	    (zero-or-more " ")
	    "="
	    (zero-or-more " ")
	    (submatch (one-or-more (not (in " ;")))) ; n
	    (zero-or-more " ")
	    ";"
	    (zero-or-more " ")
	    (one-or-more (not (in " "))) ; x
	    "<="
	    (zero-or-more " ")
	    (submatch (one-or-more (not (in " ;")))) ; m
	    (zero-or-more " ")
	    ";"
	    (zero-or-more " ")
	    (one-or-more (not (in " =+"))) ; x
	    (zero-or-more " ")
	    "+="
	    (zero-or-more " ")
	    "1"
	    (zero-or-more " ")
	    ")"
))
       "for \\1 in range(\\2, \\3 + 1):")

  ;; Change foreach (x as y) to for y in x:
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "foreach"
	    (zero-or-more " ")
	    "("
	    (zero-or-more " ")
	    (submatch (one-or-more (not-char " "))) ; x
	    " as "
	    (submatch (one-or-more (not (in " )")))) ; y
	    (zero-or-more " ")
	    ")"))
   "for \\2 in \\1:")
	    
  ;; Remove parentheses around if conditions
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "if ("
	    (submatch (minimal-match (zero-or-more anything)))
	    ")"))
   "if \\1:")

  ;; Change elseif { to elif:
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "elseif"
	(minimal-match (zero-or-more anything))
	"{"))
   "elif:")

  ;; Change else { to else:
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "else"
	(minimal-match (zero-or-more anything))
	"{"))
   "else:")

  ;; Change functions to defs
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "function "
	(submatch (minimal-match (zero-or-more anything)))
	")"
	(minimal-match (zero-or-more anything))
	"{"))
   "def \\1):")

  ;; Change class syntax
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "class "
	    (submatch (one-or-more not-newline))))
   "class \\1:")

  ;; Change class variable declarations
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "var "
	    (submatch (zero-or-more (not-char " ")))))
       "\\1 = None")


  ;; Change is_array(x) to isinstance(x, list)
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "is_array"
	    (zero-or-more (char " "))
	    (char "(")
	    (zero-or-more (char " "))
	    (submatch (zero-or-more (not-char " ")))
	    (zero-or-more (char " "))
	    (char ")")))
   "isinstance(\\1, list)")


  ;; Change arrays to lists
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "array(" 
	    (zero-or-more (char " "))
	    (submatch (zero-or-more (not-char ")")))
	    (char ")")))
   "[\\1]")

  ;; Change count(x) to len(x)
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "count(" 
	    (submatch (zero-or-more (not-char ")")))
	    (char ")")))
   "len(\\1)")

  ;; Change print_r(x, true) to str(x)
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "print_r(" 
	    (submatch (zero-or-more (not-char ")")))
	    (char ")")
	    (zero-or-more (char " "))
	    (char ",")
	    (zero-or-more (char " "))
	    "true"
	    (zero-or-more (char " "))
	    (char ")")))
   "str(\\1)")

  ;; Change print_r(x) to print str(x)
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "print_r(" 
	    (submatch (zero-or-more (not-char ")")))
	    (char ")")))
   "print str(\\1)")


  ;; Change class instance calls
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "new "
	    (submatch (one-or-more (not-char "(")))
	    "("
	    (submatch (zero-or-more (not-char ")")))
	    ")"))
       "\\1(\\2)")

  ;; Change this to self
  (goto-char (point-min))
  (query-replace-regexp 
   (rx (seq "this" (char ".")))
   "self.")

  ;; Get rid of brackets
  (goto-char (point-min))
  (query-replace-regexp (rx (or "{" "}")) "")

  ;; Change is_numeric(x) to x.is_digit()
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "is_numeric"
	    (zero-or-more (char " "))
	    (char "(")
	    (zero-or-more (char " "))
	    (submatch (zero-or-more (not-char " ")))
	    (zero-or-more (char " "))
	    (char ")")))
   "\\1.isdigit()")

  ;; Change strpos(x) to x.find()
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "strpos"
	    (zero-or-more (char " "))
	    (char "(")
	    (zero-or-more (char " "))
	    (submatch (zero-or-more (not-char " ")))
	    (zero-or-more (char " "))
	    (char ",")
	    (zero-or-more (char " "))
	    (submatch (zero-or-more (not-char " ")))
	    (zero-or-more (char " "))
	    (char ")")))
   "\\1.find(\\2)")


  ;; Change array_key_exists(key,array) to array.has_key(key)
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "array_key_exists"
	    (zero-or-more (char " "))
	    (char "(")
	    (zero-or-more (char " "))
	    (submatch (zero-or-more (not-char " ")))
	    (zero-or-more (char " "))
	    (char ",")
	    (zero-or-more (char " "))
	    (submatch (zero-or-more (not-char " ")))
	    (zero-or-more (char " "))
	    (char ")")))
   "\\2.has_key(\\1)")


  ;; Change echo to print
  (goto-char (point-min))
  (query-replace-regexp "echo" "print")

  ;; Change trigger_error to raise
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "trigger_error("
	    (submatch (zero-or-more (not-char ")")))	    
	    (char ")")))
   "raise StandardError(\\1)")

  ;; Change false to False
  (goto-char (point-min))
  (query-replace-regexp "false" "False")

  ;; Change true to True
  (goto-char (point-min))
  (query-replace-regexp "true" "True")

  ;; Change NULL to None
  (goto-char (point-min))
  (query-replace-regexp "NULL" "None")

  ;; Change None!=x to x is not None
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "None"
	    (zero-or-more (char " "))
	    "!="
	    (zero-or-more (char " "))
	    (submatch (zero-or-more (not-char " ")))))
   "\\1 is not None")

  ;; Change x!=None to x is not None
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq (zero-or-more (char " "))
	    (submatch (zero-or-more (not-char " ")))
            (zero-or-more (char " "))
	    "!="
            (zero-or-more (char " "))
	    "None"))
   "\\1 is not None")

  ;; Change None==x to x is None
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "None"
	    (zero-or-more (char " "))
	    "=="
	    (zero-or-more (char " "))
	    (submatch (zero-or-more (not-char " ")))))
   "\\1 is None")

  ;; Change x==None to x is None
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq (zero-or-more (char " "))
	    (submatch (zero-or-more (not-char " ")))
            (zero-or-more (char " "))
	    "=="
            (zero-or-more (char " "))
	    "None"))
   "\\1 is None")

  ;; Change RoboDoc todos to Epydoc todos
  (goto-char (point-min))
  (query-replace-regexp "@todo " "@todo: ")

  ;; Change RoboDoc notes to Epydoc notes
  (goto-char (point-min))
  (query-replace-regexp "@note " "@note: ")

  ;; Change RoboDoc versions to Epydoc versions
  (goto-char (point-min))
  (query-replace-regexp "@version " "@version: ")

  ;; Change RoboDoc parameters to Epydoc parameters
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "@param"
	    " "
	    (submatch (one-or-more (not-char " "))) ; type
	    " "
	    (submatch (one-or-more (not-char " "))) ; name
	    " "
	    (submatch (minimal-match (one-or-more anything)))
	    line-end))
   "@param \\2: \\3\n\t@type \\2: \\1")

  ;; Change RoboDoc return values to Epydoc return values
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "@return"
	    " "
	    (submatch (one-or-more (not-char " "))) ; type
	    " "
	    (submatch (minimal-match (one-or-more anything))) ; desc
	    line-end))
   "@return: \\2 \n\t@rtype: \\1")

  ;; Change WARNING to Epydoc warning
  (goto-char (point-min))
  (query-replace-regexp "WARNING: " "@warning: ")

  ;; Move docstrings under function definitons
  (goto-char (point-min))
  (query-replace-regexp
   (rx (seq "\"\"\""
	    (submatch (minimal-match (zero-or-more (not-char "\""))))
	    "\"\"\"\ndef"
	    (submatch (minimal-match (zero-or-more anything)))
	    ":"))
   "def \\2:\n\"\"\"\\1\"\"\"")

) ; End defun php2py

(provide 'php2py)

;; $Id: php2py.el,v 1.2 2006/11/23 12:09:06 aarre Exp $
;;; php2py.el ends here