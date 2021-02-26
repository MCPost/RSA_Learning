
% SVM and LDA
cfg = [];
cfg.metric      = 'acc';
cfg.cv          = 'leaveout';
cfg.repeat      = 5;
cfg.classifier  = 'lda';
cfg.feedback	= 0;
[lda_auc_class] = mv_classify(cfg, Cur_trial, [ones(1,size(cur_trial1,1)) 2*ones(1,size(cur_trial2,1))]');
