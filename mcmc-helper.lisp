(defpackage mcmc-helper
  (:use :cl)
  (:local-nicknames (:gsl :gsll)
                    (:plt :eazy-gnuplot))
  (:export
   ;; 一些注记用的函数
   #:square #:random-sign #:random-vect 
   ;; 一些 list 统计函数
   #:sum #:average #:average #:var #:hist #:table-hist #:mc-integrate-hist-table
   #:average-table #:variance-table #:at #:near
   ;; 一些 list 向量函数和矩阵相关函数
   #:vec-length #:vec-num-product #:vec-vec-dot-product
   #:table-transpose
   ;; 一些绘图函数
   #:plot-table #:plot-hist-table #:plot-grid))

(in-package mcmc-helper)

(defun square (x) (* x x))

(defun sum (lst) (coerce (apply #'+ lst) 'float))

(defun average (lst) (coerce (/ (sum lst) (length lst)) 'float))

(defun var (lst)
  "var(LST) = average(LST^2) - average(LST)^2"
  (- (average (mapcar #'square lst)) (square (average lst))))

(defun vec-length (vec)
  "一个向量的长度: |vec| = √(∑ vec_i^2)."
  (sqrt (apply #'+ (mapcar #'square vec))))

(defun vec-num-product (vec num)
  "向量的数乘: num x vec."
  (loop for vi in vec collect (* num vi)))

(defun vec-vec-dot-product (vec1 vec2)
  "向量的点乘: vec1 · vec2"
  (apply #'+ (mapcar #'* vec1 vec2)))

(defun table-transpose (table)
  "转置矩阵"
  (if (atom (car table))
      (loop for elem in table collect (list elem))
      (apply #'mapcar #'list table)))

(defun hist (dat &key (bin-num 100) (x-min 'auto) (x-max 'auto))
  "从列表数据中生成 `bin-num' 个数据 bin 的范围为 `x-min' 到 `x-max' 的直方统计."
  (let ((hist (gsl:make-histogram bin-num)))
    (gsl:set-ranges-uniform
     hist
     (coerce (if (eq x-min 'auto) (apply #'min dat) x-min) 'double-float)
     (coerce (if (eq x-max 'auto) (apply #'max dat) x-max) 'double-float))
    (loop for record in dat do
      (gsl:increment hist (coerce record 'double-float)))
    hist))

(defun mc-integrate-hist (func hist)
  "把函数积分转换为统计采样和平均."
  (let ((x-min    (gsl:min-range hist))
        (x-max    (gsl:max-range hist))
        (bin-num  (first (grid:dimensions hist)))
        (total    0)
        (integral 0))
    (loop for i below bin-num
          for x from x-min by (/ (- x-max x-min) bin-num)
          do (setf total (+ total (grid:aref hist i)))
          do (setf integral (+ integral (* (funcall func x) (grid:aref hist i)))))
    (/ integral total)))

(defun table-hist (table &key (heading T) (bin-num 100) (x-min 'auto) (x-max 'auto))
  "从表格数据中提取直方图表."
  (if heading
      (list (first table)
            (loop for dat in (table-transpose (rest table))
                  collect (hist dat :bin-num bin-num :x-min x-min :x-max x-max)))
      (list (loop for dat in (table-transpose table)
                  collect (hist dat :bin-num bin-num :x-min x-min :x-max x-max)))))

(defun mc-integrate-hist-table (table &key (heading T) (func #'identity))
  "对表格对应的直方图进行积分."
  (if heading
      (list (first table)
            (loop for hist in (second table)
                  collect (mc-integrate-hist func hist)))
      (loop for hist in (second table)
            collect (mc-integrate-hist func hist))))

(defun maptable (table func &key (heading T))
  "对表格元素使用 `func' 进行历遍."
  (if heading
      (list (first table)
            (loop for item in (table-transpose (rest table))
                   collect (funcall func item)))
      (loop for item in (table-transpose (rest table))
             collect (funcall func item))))

(defun average-table (table &key (heading T))
  "对表格元素求平均."
  (maptable table #'average :heading heading))

(defun variance-table (table &key (heading T))
  "对表格元素求方差."
  (maptable table #'var :heading heading))

(defmacro at (mat &rest pos)
  "方便访问结构嵌套中的元素的宏.
例: `(at model i j)' 会被展开为 `(nth j (nth i model))'."
  (if (null pos) mat `(nth ,(car pos) (at ,mat . ,(cdr pos)))))

(defmacro near (model &rest pos)
  "对于嵌套列表读取周围元素, 返回一个临近元素列表, 使用周期性边界条件.
例: `(near model (i n) (j m))' 会被展开为 `(list (at model (mod (1+ i) n) j) ...)'."
  (let ((idx-pos (mapcar #'first pos)))
    (cons
     'list
     (loop for (idx+ idx-) in (loop for (idx n) in pos
                                    collect `((mod (1+ ,idx) ,n) (mod (1- ,idx) ,n)))
           for i from 0
           for before = (subseq idx-pos 0 i)
           for after  = (subseq idx-pos (1+ i))
           collect (append `(at ,model) before (list idx+) after)
           collect (append `(at ,model) before (list idx-) after)))))

(defun random-sign (&optional (positive 0.5))
  "按照正数出现的概率为 `positive' 随机生成 +1 或者 -1."
  (if (< (random 1.0) positive) +1 -1))

(defun random--vect (dimension)
  "返回一个 `dimension' 维的长度不定的随机向量, 各分量 -1 ~ 1 取值."
  (loop for - below dimension collect (* (random-sign) (random 1.0))))

(defun random-vect (&key (dimension 2) (vec-length 1) (threshold 1e-5))
  "返回一个 `dimension' 维的长度为 `vec-length' 的随机向量.
+ `vec-length' 为 `random' 时, 结果同 `random--vect',
+ `threshold'  为防止随机生成向量太短而导致归一化时出错的保险."
  (if (eq vec-length 'random)
      (random--vect dimension)
      (let ((vec (loop for temp = (random--vect dimension)
                       while (<= (vec-length temp) threshold)
                       finally (return temp))))
        (vec-num-product vec (/ vec-length (vec-length vec))))))

(defun plot-table (output table
                   &key (heading T) (debug NIL) (font "Arial")
                     (idx-time #'identity) (plot-style :lines)
                     (multiplot NIL))
  "根据表的信息绘制图像."
  (let ((plts (table-transpose table)))
    (plt:with-plots (s :debug debug)
      (plt:gp-setup :terminal '(pngcairo) :output output :font font)
      (if multiplot (plt:gp :set :multiplot (list multiplot)))
      (loop for plt in plts do
        (plt:plot (lambda ()
                    (loop for dat in (if heading (rest plt) plt)
                          for idx from 1
                          do (format s "~&~a ~a" (funcall idx-time idx) dat)))
                  :with (if (listp plot-style) plot-style (list plot-style))
                  :title (if heading (string (first plt)) "")))
      (if multiplot (plt:gp :unset 'multiplot))))
  output)

(defun plot-hist-table (output table &key (debug NIL) (heading T))
  "根据直方图表格进行绘制."
  (let ((plts (table-transpose table)))
    (print plts)
    (plt:with-plots (s :debug debug)
      (plt:gp-setup :terminal '(pngcairo) :output output)
      (loop for histo in plts do
        (plt:plot (lambda ()
                    (let* ((hist    (if heading (second histo) histo))
                           (x-min   (gsl:min-range   hist))
                           (x-max   (gsl:max-range   hist))
                           (bin-num (first (grid:dimensions hist))))
                      (loop for idx below bin-num
                            for x from x-min by (/ (- x-max x-min) bin-num)
                            do (format s "~&~,4f ~,1f" x (grid:aref hist idx)))))
                  :with '(:boxes)
                  :title (if heading (string (first histo)) ""))))
    output))

(defun plot-grid (output grid
                       &key (debug NIL))
  "绘制二维表格/矩阵的图片."
  (plt:with-plots (s :debug debug)
    (plt:gp-setup :terminal '(pngcairo) :output output)
    (plt:gp :set :autoscale :xfix)
    (plt:gp :set :autoscale :yfix)
    (plt:gp :set :autoscale :cbfix)
    (plt:gp :unset :colorbox)
    (plt:gp :set :palette :grey)
    (plt:plot (lambda ()
                (loop for row in grid do (format s "~&~{~a~^ ~}" row)))
              :matrix '() :with '(image) :notitle '()))
  output)
