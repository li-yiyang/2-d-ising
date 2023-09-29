(defpackage metropolis
  (:use :cl :mcmc-helper)
  (:export #:sampling))

(in-package metropolis)

;;; 通用的 Metropolis 采样函数主体
(defun sampling (init update samples 
                 &key (burn-in 0) (collect #'identity) (sample-step 1))
  "采样函数, 初始状态 `init', 
   每次根据当前状态更新一个新的状态 `update' 函数, 
   采样 `samples' 个样本: 
+ 更新状态的函数 `update' 需要接受一个当前状态作为参数, 且返回值为下一个状态. 
+ 初始抛弃采样点数量 `burn-in' 个
+ 采样时每 `sample-step' 个采样点进行一次数据采集
+ 采集的数据通过 `collect' 函数处理后保存为列表返回"
  (let ((state init))
    (loop for n below burn-in do
      (setf state (funcall update state)))
    (loop for n below samples do
      (loop for s below sample-step do
        (setf state (funcall update state)))
          collect (funcall collect state))))
