(defpackage 2-d-ising
  (:use :cl :mcmc-helper :metropolis))

(in-package 2-d-ising)

;;; 定义 2 维 Ising 模型的能量, 磁矩, 每次更新的翻转函数和每次仿真的记录函数
;;; 模型描述
(defun 2-d-ising-model (n m initial-state)
  "二维 Ising 模型的数据描述."
  (cond ((functionp initial-state)
         (loop for j below m
               collect (loop for i below n
                             collect (funcall initial-state i j))))
        ((eq initial-state 'random)
         (2-d-ising-model
          n m (lambda (i j) (declare (ignore i j)) (1- (* 2 (random 2))))))
        (T (2-d-ising-model
            n m (lambda (i j) (declare (ignore i j)) initial-state)))))

;;; 模型磁矩
(defun 2-d-ising-m (model)
  "计算 2 维 Ising 模型的磁矩."
  (sum (mapcar #'sum model)))

;;; 模型能量
(defun 2-d-mag-field (field)
  "根据 `field' 描述符生成一个二维磁场."
  (cond ((eq field 'unify)     (lambda (i j) (declare (ignore i j)) 1))
        ((eq field 'neg-unify) (lambda (i j) (declare (ignore i j)) -1))
        ((eq field 'none)      (lambda (i j) (declare (ignore i j)) 0))
        ((functionp field)     field)
        (T                     (lambda (i j) (declare (ignore i j)) field))))

(defun 2-d-ising-e (model &key (field 'none))
  "二维 Ising 模型的能量计算."
  (let ((energy 0)
        (mag-field (2-d-mag-field field))
        (n (length (first model)))
        (m (length model)))
    (loop for j below m do
      (loop for i below n do
        (setf energy
              (- energy
                 ;; 互能项, 剔除重复项除以 4 
                 (/ (* (at model i j) (sum (near model (i n) (j m)))) 4)
                 ;; 外场项
                 (* (at model i j) (funcall mag-field i j))))))
    energy))

;;; 更新函数和收集函数
(defun 2-d-ising-gen-update (model n m kT &key (field 'unify) (give-up 'never))
  "生成一个 2 维 Ising 模型的更新函数. 
为了加速更新速度并且减少中间量计算开销, 使用每次翻转的位置 `(i . j)' 作为状态.
并且会修改变量 model.

对于 dE 的计算: 因为 E = - ∑ Si Sj - ∑ Hi Si, 所以在将 Sk 变为 -Sk 的时候,
导致的能量变化为: dE = 2(∑ Sk Sj + Hk Sk)."
  (let ((mag-field (2-d-mag-field field)))
    (if (eq give-up 'never)
        ;; 不考虑放弃, 容易在温度低时一直尝试生成
        (lambda (state) (declare (ignore state))
          (loop for i = (random n) for j = (random m)
                while (let* ((dE (* 2 (at model i j)
                                    (+ (sum (near model (i n) (j m)))
                                       (funcall mag-field i j))))
                             (rand (random 1.0))
                             (p (exp (- (/ dE kT)))))
                        ;; 如果能量变化小于零, 或者概率小于接受概率就终止
                        (not (or (< dE 0) (< rand p))))
                finally (progn (setf (at model i j) (- (at model i j)))
                               (return (cons i j)))))
        ;; 考虑放弃, 在 `give-up' 次内进行尝试
        (lambda (state) (declare (ignore state))
          (loop for i = (random n) for j = (random m)
                for count from 0
                while (let* ((dE (* 2 (at model i j)
                                    (+ (sum (near model (i n) (j m)))
                                       (funcall mag-field i j))))
                             (rand (random 1.0))
                             (p (exp (- (/ dE kT)))))
                        ;; 如果能量变化小于零, 或者概率小于接受概率就终止
                        (and (not (or (< dE 0) (< rand p)))
                             (< count give-up)))
                finally (progn (setf (at model i j) (- (at model i j)))
                               (return (cons i j))))))))
(defun 2-d-ising--gen-collector (model collect-fn)
  "生成一个收集函数, `collect-fn' 需为一个函数列表."
  (lambda (state) (declare (ignore state))
    (loop for fn in collect-fn
          collect (funcall fn model))))

(defun 2-d-ising-gen-collector (model collect &key (field 'unify))
  "生成一个收集函数."
  (let ((collector-fn `((energy . ,(lambda (state)
                                     (declare (ignore state))
                                     (2-d-ising-e model :field field)))
                        (mag-m  . ,#'2-d-ising-m)
                        (state  . ,#'identity)
                        (model  . ,#'(lambda (state)
                                       (declare (ignore state)) (copy-tree model))))))
    (if (atom collect)
        (2-d-ising-gen-collector model (list collect))
        (2-d-ising--gen-collector
         model (mapcar (lambda (item) (cdr (assoc item collector-fn))) collect)))))

;;; 对于 2 维 Ising 模型的计算函数
(defun 2-d-ising-sim (n samples
                      &key (m n) (initial-state 1) (kT 1) (field 'unify)
                        (sample-step 1) (burn-in 0) (collect 'state) (give-up 'never))
  "对 2 维 Ising 模型进行一个真的仿."
  (let* ((model     (2-d-ising-model n m initial-state))
         (update    (2-d-ising-gen-update model n m kT :field field :give-up give-up))
         (collector (2-d-ising-gen-collector model collect)))
    (cons (if (listp collect) collect (list collect))
          (sampling model update samples
                    :sample-step sample-step :burn-in burn-in
                    :collect collector))))

;;; 一些数据收集的函数
(defun 2-d-ising-sim-collect-average
    (n kT collect &key (m n) (initial-state 1) (field 'unify) (sample-step 10)
                    (burn-in 1000) (samples 1000) (heading T) (give-up 'never))
  "将 2 维 Ising 模型的结果求平均."
  (let ((table (average-table
                (2-d-ising-sim n samples
                               :m m :initial-state initial-state
                               :kT kT :field field :sample-step sample-step
                               :burn-in burn-in :collect collect :give-up give-up))))
    (if heading table (rest table))))

(defun 2-d-ising-sim-collect-var
  (n kT collect &key (m n) (initial-state 1) (field 'unify) (sample-step 10)
                    (burn-in 1000) (samples 1000) (heading T) (give-up 'never))
  "将 2 维 Ising 模型结果求方差."
  (let ((table (variance-table
                (2-d-ising-sim n samples
                               :m m :initial-state initial-state
                               :kT kT :field field :sample-step sample-step
                               :burn-in burn-in :collect collect :give-up give-up))))
    (if heading table (rest table))))
(defun 2-d-ising-sim-collect-heat-capacity
    (n kT &key (m n) (initial-state 1) (field 'unify) (sample-step 10)
            (burn-in 1000) (samples 1000) (give-up 'never))
  "计算热容."
  (let* ((sim (2-d-ising-sim n samples :m m :field field :sample-step sample-step
                                       :kT kT :initial-state initial-state :give-up give-up
                                       :collect '(energy) :burn-in burn-in))
         (var-E (car (second (variance-table sim)))))
    (/ var-E (square kT) (square (* n m)))))
