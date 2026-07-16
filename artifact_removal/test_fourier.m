c1_path = '/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/preprocessed_CAP_data/20221201_2/CAPData_ch1.mat';
c2_path = '/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/preprocessed_CAP_data/20221201_3/CAPData_ch1.mat';
c1_data = load(c1_path);
c2_data = load(c2_path);
[~,AMP_IDX_1] = max(c1_data.amp_values_actual);
c1 = c1_data.ch_neural_avg(:, AMP_IDX_1);
[~,AMP_IDX_2] = max(c2_data.amp_values_actual);
c2 = c2_data.ch_neural_avg(:, AMP_IDX_2);

dt = c1_data.t_window(2) - c1_data.t_window(1); 

remove_artifact_fourier(c1,c2,1/dt);