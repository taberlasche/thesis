%%
%% Vorschlag für eine LaTeX2e-Dokumentenklasse für Abschlussarbeiten am Physik-Department der TU München
%% Stand 24.04.2018 (PHvL)
%%
\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{tumphthesis}[2018/08/24 LaTeX2e TUM PH Thesis layout]
% Basisklasse: KOMA-Buch (scrbook)
\PassOptionsToPackage{headinclude=true,footinclude=false}{typearea}
\LoadClass[version=last,%
	headinclude=true%
	,footinclude=false%
	,fontsize=11pt%
	,twocolumn=false%
	,twoside=false%
	,openany%
	,footnotes=multiple%
	,chapterprefix=true%
	,appendixprefix=true%
	,numbers=noenddot%
	,bibliography=totocnumbered%
	]{scrbook}%[2010/06/17]
\DefineFamily{TUMPH}
\DefineFamilyMember{TUMPH}
\def\tumph@font{9}
\FamilyNumericalKey{TUMPH}{font}{tumph@font}{%
	{cm}{9},{latex}{9},%
	{palatino}{0},{mathpazo}{0},%
	{cmsuper}{1},%
	{minionpro}{2},%
	{helvetica}{3},{helvet}{3},%
	{palatino-tum}{4},{mathpazo-tum}{4},%
	{tum}{5},%
	{lmodern}{6}%
	}
\FamilyBoolKey{TUMPH}{draft}{@draftmode}
\FamilyBoolKey{TUMPH}{dissertation}{@dissertation}
\FamilyBoolKey{TUMPH}{pstricks}{@loadpst}
\FamilyBoolKey{TUMPH}{bibunits}{@loadbibunits}
\FamilyBoolKey{TUMPH}{siunitx}{@loadsiunitx}
\FamilyBoolKey{TUMPH}{nohyperref}{@noloadhyperref}
\FamilyBoolKey{TUMPH}{nopdfcomment}{@noloadpdfcomment}
\FamilyBoolKey{TUMPH}{nocsquotes}{@noloadcsquotes}
\FamilyBoolKey{TUMPH}{addfonts}{@loadadditionalfonts}
\FamilyBoolKey{TUMPH}{theorem}{@loadtheorem}
\newcommand{\file@encoding}{utf8}
\DefineFamilyKey{TUMPH}{encoding}{\renewcommand{\file@encoding}{#1}}
\FamilyProcessOptions{TUMPH}\relax
\if@draftmode
	\relax
\else
	\KOMAoptions{BCOR=1cm%
		,DIV=10%
		,open=any%
		,fontsize=11pt%
		,twoside=true%
		}
\fi
\if@noloadhyperref
	\FamilySetBool{TUMPH}{nopdfcomment}{@noloadpdfcomment}{true}
\else
\PassOptionsToPackage{\file@encoding}{inputenc}
\RequirePackage{inputenc}
\PassOptionsToPackage{T1}{fontenc}
\RequirePackage{fontenc}
\PassOptionsToPackage{ngerman,british}{babel}
\RequirePackage{babel}
\selectlanguage{ngerman}
\if@noloadcsquotes
	\relax
\else
	\PassOptionsToPackage{german=guillemets}{csquotes}
	\RequirePackage{csquotes}
	\defineshorthand{"`}{\openautoquote}
	\defineshorthand{"'}{\closeautoquote}
\fi

% xcolor dient zur vereinfachten Definition von Farben.
\RequirePackage{xcolor}

% Ist die Einbindung von Postscript Grafiken gewünscht, so geschieht dies
% durch pst-pdf.
\if@loadpst
	\PassOptionsToPackage{nosub}{psfragx}
	\RequirePackage{psfragx}
	\PassOptionsToPackage{notightpage}{pst-pdf}
	\RequirePackage{pst-pdf}
	\if@draftmode
		\setkeys{Gin}{showname=true}
	\fi
\else
	\RequirePackage[final]{graphicx}
\fi

% Um PDF-Kommentare für ToDos o. Ae. zu setzen, wird pdfcomment verwendet
% und ein \todo-Makro definiert.
\if@noloadpdfcomment
	\newcommand{\todo}[1]{%
		\marginpar{#1}%
		\if@draftmode
			\relax
		\else
			\ClassWarning{tumphthesis}{Warning! There are \string\todo\space notes in the final document.}
		\fi}
\else
	\RequirePackage{pdfcomment}
	\newcommand{\todo}[1]{%
		\pdfcomment[color=red,open=true]{#1}%
		\if@draftmode
			\relax
		\else
			\ClassWarning{tumphthesis}{Warning! There are \string\todo\space notes in the final document.}
		\fi}
\fi

% biblatex wird zum Setzen des Literatirverzeichnisses verwendet
\PassOptionsToPackage{sorting=none,autolang=other,backend=biber,
	bibstyle=numeric-comp,citestyle=numeric-comp,%
	}{biblatex}
\RequirePackage{biblatex}
\IfFileExists{\jobname.bib}{\addbibresource{\jobname.bib}}{\relax}%
\defbibheading{bibliography}[\bibname]{\chapter{#1}}

% Die AMS-Pakete dienen zum modernen Satz von Gleichungen. Die Integralgrenzen
% werden grundsätzlich unter/über das Integralzeichen gesetzt.
\PassOptionsToPackage{intlimits}{amsmath}
\RequirePackage{amsmath}
\RequirePackage{amstext}
\RequirePackage{nicefrac}
\PassOptionsToPackage{showonlyrefs}{mathtools}
\RequirePackage{mathtools}
\RequirePackage{amsfonts}
\RequirePackage{newclude}
\RequirePackage{bm}
\RequirePackage{listings}
% Sätze, Definitionen und Beweise können mit amsthm und shadethm einfach und
% sauber gesetzt werden.
\if@loadtheorem
	\RequirePackage{amsthm}
	\RequirePackage{shadethm}
	\newshadetheorem{hauptsatz}{Hauptsatz}
	\newshadetheorem{satz}{Satz}
	\newshadetheorem{hilfssatz}{Hilfssatz}
	\newshadetheorem{theorem}{Theorem}
	\newshadetheorem{lemma}{Lemma}
	\newshadetheorem{folgerung}{Folgerung}
	\newshadetheorem{proposition}{Proposition}
	\theoremstyle{definition}
	\newshadetheorem{definition}{Definition}
	\newshadetheorem{beispiel}{Beispiel}
	\newshadetheorem{bemerkung}{Bemerkung}
	\newshadetheorem{axiom}{Axiom}
	\AtBeginDocument{\setlength{\shadedtextwidth}{\linewidth}}
\fi

\RequirePackage{booktabs}
\RequirePackage{caption}
\RequirePackage{multirow}
\RequirePackage{xspace}
\RequirePackage{float}
% Mit dem Paket siunitx lassen sich einheitenbehaftete Größen auf einfache Weise
% korrekt setzen.
\if@loadsiunitx
	\PassOptionsToPackage{locale=DE,per-mode=symbol-or-fraction,load-configurations={binary,abbreviations}}{siunitx}
	\RequirePackage{siunitx}[2010/07/09]
\fi

% Ein paar Schriften für Sonderzeichen.
\if@loadadditionalfonts
	\RequirePackage{amsfonts}
	\RequirePackage{mathrsfs}
	\RequirePackage{textcomp}
	\RequirePackage{dsfont}
	\newcommand{\IC}{\ensuremath{\mathds{C}}}
	\newcommand{\IK}{\ensuremath{\mathds{K}}}
	\newcommand{\IN}{\ensuremath{\mathds{N}}}
	\newcommand{\IQ}{\ensuremath{\mathds{Q}}}
	\newcommand{\IR}{\ensuremath{\mathds{R}}}
	\newcommand{\IZ}{\ensuremath{\mathds{Z}}}
	\newcommand{\1}{\ensuremath{\mathds{1}}}
\else
	\newcommand{\IC}{\ensuremath{\mathbb{C}}}
	\newcommand{\IK}{\ensuremath{\mathbb{K}}}
	\newcommand{\IN}{\ensuremath{\mathbb{N}}}
	\newcommand{\IQ}{\ensuremath{\mathbb{Q}}}
	\newcommand{\IR}{\ensuremath{\mathbb{R}}}
	\newcommand{\IZ}{\ensuremath{\mathbb{Z}}}
	\newcommand{\1}{\ensuremath{\mathbb{I}}}
\fi

% Einstellung der Brotschrift
\def\Enable@MicroType{%
	\if@draftmode
		\relax
	\else
		\PassOptionsToPackage{babel}{microtype}
		\RequirePackage{microtype}
		\DisableLigatures{encoding=T1,family=tt*}
	\fi
}
\ifcase\tumph@font % Palatino, Pazo (Mathe), Helvetica (Sans Serif), Courier (TT) mit Mikrotypographie
	\RequirePackage{mathpazo}
	\RequirePackage{courier}
	\RequirePackage[scaled=.95]{helvet}
	\Enable@MicroType
\or % LaTeX-Standardschriften mit Mikrotypographie
	\Enable@MicroType
\or % MinionPro, MnSymbol (Mathe) mit Mikrotypographie
	\RequirePackage{MinionPro}
	\RequirePackage{MnSymbol}
	\Enable@MicroType
\or % Helvetica, Euler-VM (Mathe), Courier (TT) mit Mikrotypographie
	\RequirePackage{eulervm}
	\RequirePackage{courier}
	\RequirePackage[scaled=.95]{helvet}
	\renewcommand{\familydefault}{\sfdefault}
	\Enable@MicroType
\or % Palatino, Pazo (Mathe), TUM-Neue-Helvetica (Sans Serif), Courier (TT) mit Mikrotypographie
	\RequirePackage{mathpazo}
	\RequirePackage{courier}
	\renewcommand{\sfdefault}{lhv}
	\Enable@MicroType
\or % TUM-Neue-Helvetica, Euler-VM (Mathe), Courier (TT) mit Mikrotypographie
	\RequirePackage{eulervm}
	\RequirePackage{courier}
	\renewcommand{\sfdefault}{lhv}
	\renewcommand{\familydefault}{\sfdefault}
	%\Enable@MicroType
\or % Latin-Modern mit Mikrotypographie
	\RequirePackage{lmodern}
	\Enable@MicroType
\else % LaTeX-Standardschriften
	\relax
\fi
\if@draftmode
	\relax
\else
	\RequirePackage{tocbasic}[2010/06/09]
	\AtEndOfClass{\RequirePackage{tocstyle}%
		\usetocstyle{allwithdot}}
\fi

% Hyperref zum setzen von Links. URLs können bei jedem Zeichen umbrochen werden.
\PassOptionsToPackage{hyphens}{url}
\if@noloadhyperref
	\RequirePackage{url}
\else
	\RequirePackage{hyperref}
\fi
\g@addto@macro\UrlBreaks{
  \do\a\do\b\do\c\do\d\do\e\do\f\do\g\do\h\do\i\do\j
  \do\k\do\l\do\m\do\n\do\o\do\p\do\q\do\r\do\s\do\t
  \do\u\do\v\do\w\do\x\do\y\do\z\do\&\do\1\do\2\do\3
  \do\4\do\5\do\6\do\7\do\8\do\9\do\0}
\g@addto@macro\UrlSpecials{\do\/{\mbox{\UrlFont/}\hskip 0pt plus 1pt}}

% Abkürzungen
\newcommand{\zB}{z.\,B.\xspace}
\newcommand{\eg}{e.\,g.\xspace}
%\newcommand{\dh}{d.\,h.\xspace}

% Einige Sonderzeichen
\renewcommand{\d}{\mathrm{d}}											% d für Differentiale
\newcommand{\e}{\mathrm{e}}												% die Eulersche Zahl e
\newcommand{\iu}{\mathrm{i}}											% die imaginäre Einheit
\newcommand{\laplace}{\ensuremath{\Delta}}								% der Laplace-Operator
\newcommand{\dalembert}{\ensuremath{\Box}} 						 		% der D'Alembert-Operator
\newcommand{\speedoflight}{\ensuremath{c}}								% Symbol für die Lichtgeschwindigkeit
\newcommand{\GNewton}{G_{\text{N}}}										% Symbol für Newtons Gravitationskonstante
\newcommand{\kb}{k_{\text{B}}}											% Symbol für die Boltzmannkonstante

% Ein paar Operatoren
\DeclareMathOperator{\diag}{diag}										% Diagonalmatrix
\renewcommand{\Re}{\operatorname{Re}}                           		% Realteil
\renewcommand{\Im}{\operatorname{Im}}                    				% Imaginärteil
\newcommand{\pri}[1]{{{#1}^{\prime}}}
\newcommand{\ppri}[1]{{{#1}^{\prime\prime}}}
\newcommand{\pdot}[1]{{\dot{#1}^{\prime}}}

% Totale und partielle Ableitungen
\newcommand{\tdiff}[3][\empty]{\ifx\empty#1
		\frac{\d\,#2}{\d #3}
	\else
		\frac{\d^{#1}\,#2}{\d #3^{#1}}
	\fi}		% (#1-te) totale Ableitung von #2 nach #3
\newcommand{\pdiff}[3][\empty]{\ifx\empty#1
		\frac{\partial\,#2}{\partial #3}
	\else
		\frac{\partial^{#1}\,#2}{\partial #3^{#1}}
	\fi}		% (#1-te) partielle Ableitung von #2 nach #3
\newcommand{\vdiff}[3][\empty]{\ifx\empty#1
		\frac{\delta\,#2}{\delta #3}
	\else
		\frac{\delta^{#1}\,#2}{\delta #3^{#1}}
	\fi}		% (#1-te) Variationsableitung von #2 nach #3

% Weitere Mathe-Makros
\DeclarePairedDelimiter{\abs}{\lvert}{\rvert}								% Betrag von #1
\DeclarePairedDelimiter{\norm}{\lVert}{\rVert}								% Norm von #1
\DeclarePairedDelimiter{\gaussbracket}{\lfloor}{\rfloor}					% Gaussklammer von #1
\DeclarePairedDelimiter{\normalorder}{\mathbf{:}}{\mathbf{:}}				% Normalordnung von #1
\newcommand{\evalat}[2]{\ensuremath{\left. {#1} \right|_{#2}}}				% #1 ausgewertet an der Stelle #2
\newcommand{\cc}[1]{\ensuremath{\overline{#1}}}								% das komplex Konjugierte von #1
\newcommand{\melement}[3]{\left<#1\left|\vphantom{#1}#2\vphantom{#3}\right|#3\right>}	% das Matrixelement <#1|#2|#3>
\newcommand{\braket}[2]{\left<\left.#1\vphantom{#2}\right|#2\right>}		% das Skalarprodukt <#1|#2>
\newcommand{\bra}[1]{\left<#1\right|}										% der Bra-Vektor <#1|
\newcommand{\ket}[1]{\left|#1\right>}										% der Ket-Vektor |#1>
\newcommand{\vev}[1]{\ensuremath{\melement{0}{#1}{0}}}						% der Vakuumerwartungswert <0|#1|0>
\newcommand{\expvalue}[2]{\ensuremath{\melement{#2}{#1}{#2}}}				% der Erwartungswert <#2|#1|#2>
\newcommand{\poissonbracket}[2]{\ensuremath{\left\{#1;#2\right\}_{\text{P.B.}}}}	% die Poissonklammer {#1;#2}
\newcommand{\antipoissonbracket}[2]{\ensuremath{\left\{#1;#2\right\}_{\text{P.B.},+}}}	% die Anti-Poissonklammer {#1;#2}_+
\newcommand{\commutator}[2]{\ensuremath{\left[#1;#2\right]}}				% der Kommutator von #1 und #2
\newcommand{\anticommutator}[2]{\ensuremath{\left\{#1;#2\right\}}}			% der Anti-Kommutator von #1 und #2
\newcommand{\susycommutator}[2]{\ensuremath{\left[#1;#2\right\}}}			% der SuSy-Kommutator von #1 und #2
\newcommand\varvec{} \let\varvec=\vec										% das Vektorsymbol mit Pfeil ist über \varvec weiter verfügbar
\renewcommand{\vec}[1]{\ensuremath{\boldsymbol{#1}}}						% Vektoren werden durch Fettdruck gekennzeichnet
\newcommand{\oplusint}{\oplus\mkern-19mu\int}								% ein Integral mit überlagertem \oplus
\newcommand{\fmslash}[2][0mu]{%
	\mathchoice
		{\fmsl@sh\displaystyle{#1}{#2}}%
		{\fmsl@sh\textstyle{#1}{#2}}%
		{\fmsl@sh\scriptstyle{#1}{#2}}%
		{\fmsl@sh\scriptscriptstyle{#1}{#2}}}
\newcommand{\fmsl@sh}[3]{%
	\m@th\ooalign{$\hfil#1\mkern#2/\hfil$\crcr$#1#3$}}						% der Feynman-Slash

\if@noloadhyperref
	\relax
\else
	\AtBeginDocument{\hypersetup{colorlinks=true,urlcolor=blue,pdfpagelayout=TwoColumnRight,pdfauthor={\@author},%
		pdftitle={\@title},%
		pdfsubject={\@subject},%
		pdfkeywords={Technische Universit\"at M\"unchen;Fakult\"at f\"ur Physik;\@subject}}}
\fi
\renewcaptionname{ngerman}{\proofname}{Beweis}

% Kolumnentitel im Kopf mit Linie werden mit scrpage2 eingerichtet.
\if@draftmode
	\relax
\else
	\RequirePackage[ilines,headsepline]{scrlayer-scrpage}
	\pagestyle{scrheadings}
	\newcommand\OrigBackmatter{} \let\OrigBackmatter=\backmatter
	\renewcommand{\backmatter}{\OrigBackmatter \pagestyle{scrplain}}
\fi
\if@dissertation
	\RequirePackage[absolute]{textpos}
	\RequirePackage{xparse}
	\NewDocumentCommand \DepartmentOrIRC {m}{\def\@departmentorirc{#1}}
	\DepartmentOrIRC{Fakultät für Physik}
	\NewDocumentCommand \DrDegree {m}{\def\@drdegree{#1}}
	\DrDegree{der Naturwissenschaften}
	\def\@dissertationfinal{\BooleanFalse}
	\def\@erstprueferinset{\BooleanFalse}
	\def\@erweitertekommission{\BooleanFalse}
	\def\@erstprueferin{noch nicht festgelegt}
	\NewDocumentCommand \ErstPrueferin {m}{\def\@erstprueferin{#1}
		\def\@erstprueferinset{\BooleanTrue}}
	\def\@vorsitz{noch nicht festgelegt}
	\NewDocumentCommand \Vorsitz {m}{\def\@vorsitz{#1}
		\def\@dissertationfinal{\BooleanTrue}}
	\NewDocumentCommand \ZweitPrueferin {m}{\def\@zweitprueferin{#1}}
	\ZweitPrueferin{noch nicht festgelegt}
	\NewDocumentCommand \DrittPrueferin {m}{\def\@drittprueferin{#1}
		\def\@erweitertekommission{\BooleanTrue}}
\fi
\renewcommand*{\subject@font}{\usekomafont{subtitle}}%
\newcommand*{\othertitle@font}{\usekomafont{subtitle}\mdseries}%
\newcommand{\cooperators}{} \let\cooperators=\publishers%
\renewcommand*{\maketitle}{%
	\if@dissertation
	\begin{titlepage}
		\begin{textblock*}{2.87cm}(\dimexpr(\paperwidth-3.87cm),1cm)%
			\raggedleft\ifpdfoutput{\includegraphics[width=1.87cm]{tumlogo-full}}{}
		\end{textblock*}%
		{\usekomafont{subtitle}\noindent\@departmentorirc\par}%
		\vfill%
		{\titlefont\huge\noindent\@title\par}%
		\vskip 1em
		{\ifx\@subtitle\@empty\else\usekomafont{subtitle}\noindent\@subtitle\par\fi}%
		\vskip 2em
		{\othertitle@font% \lineskip 0.75em
			\noindent%\begin{tabular}[t]{c}
				\@author
			%\end{tabular}\par
		}%
		\vfill%
		{\IfBooleanTF{\@dissertationfinal}{%
			\noindent Vollständiger Abdruck der von der promotionsführenden Einrichtung
			\@departmentorirc\ der Technischen Universität München zur Erlangung
			des akademischen Grades eines Doktors \@drdegree\ genehmigten
			Dissertation.\par
			\vskip 4em
			\noindent Vorsitzende(r): \@vorsitz\par
			\vskip 1em
			\noindent Prüfende der Dissertation
			\begin{enumerate}
				\item \@erstprueferin
				\item \@zweitprueferin
				{\IfBooleanT{\@erweitertekommission}{\item \@drittprueferin}}%
			\end{enumerate}
			\vskip 2em
			Die Dissertation wurde am \@date\ bei der Technischen Universität
			München eingereicht und durch die promotionsführenden Einrichtung
			\@departmentorirc\ angenommen.}{%
			\noindent Dissertation zur Erlangung des akademischen Grades eines Doktors
			\@drdegree\ zur Vorlage bei der promotionsführenden Einrichtung
			\@departmentorirc\ der Technischen Universität München.%
			{\IfBooleanT{\@erstprueferinset}{\par\vskip 1em
			\noindent Die Dissertation wurde von \@erstprueferin\
			angeregt und in wesentlichen Teilen betreut.}}}}%
	\end{titlepage}
	\else
	\begin{titlepage}
		\noindent\begin{minipage}[b]{.5\linewidth}%
			\normalfont\sffamily Technische Universit\"{a}t M\"{u}nchen\\
			Fakult\"{a}t f\"{u}r Physik
		\end{minipage}
		\hfill
		\begin{minipage}[b]{.4\linewidth}
			\raggedleft\ifpdfoutput{\includegraphics[width=1cm]{PH}\hspace{.5cm}\includegraphics[width=1.87cm]{tumlogo}}{}
		\end{minipage}%
		\null\vfill%
		\begin{center}
			\ifx\@subject\@empty \else
				{\subject@font \@subject \par}%
				\vskip 3em
			\fi
			{\titlefont\huge \@title\par}%
			\vskip 1em
			{\ifx\@subtitle\@empty\else\usekomafont{subtitle}\@subtitle\par\fi}%
			\vskip 2em
			{\othertitle@font \lineskip 0.75em
				\begin{tabular}[t]{c}
					\@author
				\end{tabular}\par
			}%
			\vskip 1.5em
			{\othertitle@font \@date \par}%
			\vskip \z@ \@plus3fill
			{\othertitle@font \@publishers \par}%
			\vskip 3em
		\end{center}\par
		\@thanks
		\vfill\null
		\if@twoside\next@tpage
		\begin{minipage}[t]{\textwidth}
			\@uppertitleback
		\end{minipage}\par
		\vfill
		\begin{minipage}[b]{\textwidth}
			\@lowertitleback
		\end{minipage}
		\fi
		\ifx\@dedication\@empty \else
			\next@tpage\null\vfill
			{\centering \Large \@dedication \par}%
			\vskip \z@ \@plus3fill
			\if@twoside \next@tpage\cleardoubleemptypage \fi
		\fi
	\end{titlepage}
	\fi
}
\setcounter{tocdepth}{1}
\setcounter{secnumdepth}{3}
\setcapindent{0pt}

\renewcommand{\bottomfraction}{0.5}
\renewcommand{\floatpagefraction}{0.7}
\renewcommand{\dblfloatpagefraction}{0.7}

\tolerance 1414
\hbadness 1414
\emergencystretch 1.5em
\hfuzz 0.3pt
\widowpenalty=10000
\vfuzz \hfuzz
\raggedbottom
