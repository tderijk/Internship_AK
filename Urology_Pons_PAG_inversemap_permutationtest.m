"""
Created by Job van den Hurk (Scannexus, Maastricht, NL)
"""

clc

clear
close all

%% variables to change

subs = {'1008','1013','2002','2004','2011','2012','2013','2017'};
%subs = {'2002'};
rootdir = '/Volumes/research/FHML_MHeNs/Urology_Div3/Thijs de rijk/OAB fMRI/Data/';
file_prefix = 'fMRI_OAB_';
folder_structure_empty = '/fmr_files/fmr_emptybladder/';
folder_structure_full = '/fmr_files/fmr_fullbladder/';
Pons_MNI_maskfile = '/Volumes/research/FHML_MHeNs/Urology_Div3/Thijs de rijk/OAB fMRI/Data/MASK_PONS.msk';
PAG_MNI_maskfile = '/Volumes/research/FHML_MHeNs/Urology_Div3/Thijs de rijk/OAB fMRI/Data/PAG_MNI_mask_medium_v2.msk';
nrofpermutations = 1000;

% load Pons mask and compute dimensions
ponsmask = xff(Pons_MNI_maskfile);
pagmask = xff(PAG_MNI_maskfile);
for ss = 1:numel(subs)
    CSF_connectivity_results = '/Volumes/research/FHML_MHeNs/Urology_Div3/Thijs de rijk/OAB fMRI/Data/fMRI_OAB_2002/fmr_files/fmr_fullbladder/2002_resultsmap_CSFcluster_PONS_empty_full.mat';
    load(CSF_connectivity_results);
    CSF_results = abs(sort(resultsmap(:)));
    CI_99 = CSF_results(round(0.99*numel(CSF_results)));
    CI_95 = CSF_results(round(0.95*numel(CSF_results)));
    clear resultsmap
    
    disp(['Processing subject ' char(subs{ss}) '...']);
    for bb = 1:2 % empty and full bladder
        if bb == 1
            disp('Analyzing empty bladder scan...');
            subdir = [rootdir file_prefix char(subs{ss}) folder_structure_empty];
        else
            disp('Analyzing full bladder scan...');
            subdir = [rootdir file_prefix char(subs{ss}) folder_structure_full];
        end
        [files, filenames] = ListFiletypes(subdir,'_MNI.vtc',1);
        
        if numel(files) > 1
            warning(['More than one vtc file found in folder! Currently using ' char(filenames{1}) '!']);
        end
        
        
        if bb == 1
            results_prefix = 'empty';
            empty_vtc = xff(files{1});
        else
            results_prefix = 'full';
            full_vtc = xff(files{1});
        end
        vtc = xff(files{1});
        % find cluster vmp file
        
        [vmpfiles, vmpfilenames] = ListFiletypes([rootdir file_prefix char(subs{ss}) folder_structure_full],'_Clusters_MNI_medium.vmp',1);
        
        if numel(vmpfiles) > 1
            warning(['More than one vmp file found in folder! Currently using ' char(vmpfilenames{1}) '!']);
        end
        vmp = xff(vmpfiles{1});
        
        
        voxel_ix = find(ponsmask.Mask == 1);
        if bb == 1
            resultsmap = zeros([3 vmp.NrOfMaps size(ponsmask.Mask)]);
        end
        results_vmp = VMPfromVTC(empty_vtc);
        
        for cc = 1:vmp.NrOfMaps
            source_TC = mean(vtc.VTCData(:,vmp.Map(cc).VMPData>0),2);
            
            for vv = 1:numel(voxel_ix)
                target_TC = vtc.VTCData(:,voxel_ix(vv));
                R = corrcoef(source_TC,target_TC);
                
                resultsmap(bb,cc,voxel_ix(vv)) = R(2);
            end
            if cc>1
                results_vmp.Map(cc) = results_vmp.Map(1);
            end
            
            
            results_vmp.Map(cc).LowerThreshold = CI_95;
            results_vmp.Map(cc).UpperThreshold = 0.5;
            results_vmp.Map(cc).VMPData = squeeze(resultsmap(bb,cc,:,:,:));
            results_vmp.Map(cc).Name = ['Cluster ' num2str(cc)];
        end
        
        
        %         disp(['Saving results to ' char(subs{ss}) '_' results_prefix '_bladder_PAG-PONS.vmp...']);
        %         results_vmp.SaveAs([subdir char(subs{ss}) '_' results_prefix '_bladder_PAG-PONS.vmp']);
        
        
    end
    resultsmap(3,:,:,:,:) = resultsmap(2,:,:,:,:) - resultsmap(1,:,:,:,:);
    resultsmap(3,(abs(resultsmap(2,:,:,:,:))<CI_95)) = 0;
    results_vmp = VMPfromVTC(vtc);
    
    for cc = 1:vmp.NrOfMaps
        if cc>1
            results_vmp.Map(cc) = results_vmp.Map(1);
        end
        results_vmp.Map(cc).LowerThreshold = CI_95;
        results_vmp.Map(cc).UpperThreshold = 0.5;
        results_vmp.Map(cc).VMPData = squeeze(resultsmap(3,cc,:,:,:));
        results_vmp.Map(cc).Name = ['Cluster ' num2str(cc)];
    end
    %     disp(['Saving results to ' char(subs{ss}) '_full-empty_PAG-PONS.vmp...']);
    %     results_vmp.SaveAs([subdir char(subs{ss}) '_full-empty_PAG-PONS.vmp']);
    % permutations: grow random clusters and compute difference map between
    % empty and full bladder states
    modulesize = zeros(1,vmp.NrOfMaps);
    for mm = 1:vmp.NrOfMaps
        modulesize(mm) = sum(vmp.Map(mm).VMPData(:));
    end
    
    
    %if ~exist('permresultsmap','var')
    permresultsmap = zeros([nrofpermutations 3 vmp.NrOfMaps size(ponsmask.Mask)]);
    for pp = 1:nrofpermutations
        disp(['Permutation ' num2str(pp) ' / ' num2str(nrofpermutations)]);
        [clusterspace,emb] = reducevolumedims(pagmask.Mask);
        clusterspace = growmodules(size(clusterspace),clusterspace,vmp.NrOfMaps,modulesize,0);
        permclusters = zeros(size(pagmask.Mask));
        permclusters(emb(1,1):emb(1,2),emb(2,1):emb(2,2),emb(3,1):emb(3,2)) = clusterspace;
        
        for cc = 1:vmp.NrOfMaps
            for bb = 1:2
                if bb == 1
                    source_TC = mean(empty_vtc.VTCData(:,find(permclusters==cc)),2);
                else
                    source_TC = mean(full_vtc.VTCData(:,find(permclusters==cc)),2);
                end
                
                for vv = 1:numel(voxel_ix)
                    if bb == 1
                        target_TC = empty_vtc.VTCData(:,voxel_ix(vv));
                    else
                        target_TC = full_vtc.VTCData(:,voxel_ix(vv));
                    end
                    R = corrcoef(source_TC,target_TC);
                    
                    permresultsmap(pp,bb,cc,voxel_ix(vv)) = R(2);
                end
                %                 if cc>1
                %                     results_vmp.Map(cc) = results_vmp.Map(1);
                %                 end
                %
                %
                %                 results_vmp.Map(cc).LowerThreshold = CI_95;
                %                 results_vmp.Map(cc).UpperThreshold = 0.5;
                %                 results_vmp.Map(cc).VMPData = squeeze(resultsmap(bb,cc,:,:,:));
                %                 results_vmp.Map(cc).Name = ['Cluster ' num2str(cc)];
            end
            permresultsmap(pp,3,cc,:) = permresultsmap(pp,2,cc,:) - permresultsmap(pp,1,cc,:);
        end
    end
    
    save([rootdir file_prefix char(subs{ss}) '/permresultsmap_' char(subs{ss}) '.mat'],'permresultsmap', '-v7.3')
    
    
    empty_vtc.ClearObject;
    full_vtc.ClearObject;
    vmp.ClearObject;
    results_vmp.ClearObject;
    
    
    % resultsmap(3,:,:,:,:) = [];
    threshold = CI_95;
    
    voxel_ix = [];
    for ff = 1:2
        
        for cc = 1:size(resultsmap,2)
            voxel_ix = [voxel_ix; find(abs(resultsmap(ff,cc,ponsmask.Mask==1))>threshold)];
        end
        
    end
    u_voxel_ix = unique(voxel_ix);
    mask_ix = find(ponsmask.Mask == 1);
    
    % permutation assessment
    perm_correlations = zeros(size(resultsmap,2),nrofpermutations);
    for pp = 1:nrofpermutations
        ccc = 1;
        for c = 1:size(resultsmap,2)
            for cc = c+1 : size(resultsmap,2)
                r = corrcoef(squeeze(permresultsmap(pp,3,c,mask_ix(u_voxel_ix))),squeeze(permresultsmap(pp,3,cc,mask_ix(u_voxel_ix))));
                perm_correlations(ccc,pp) = r(2);
                ccc = ccc + 1;
            end
        end
    end
    
    min_perm_correlation = min(perm_correlations,[],1);
    median_perm_correlation = median(perm_correlations,1);
    mean_perm_correlation = mean(perm_correlations,1);
    
    for c = 1:size(resultsmap,2)
        for cc = c+1 : size(resultsmap,2)
            r = corrcoef(resultsmap(3,c,mask_ix(u_voxel_ix)),resultsmap(3,cc,mask_ix(u_voxel_ix)));
            p_min = ((sum(r(2)>=min_perm_correlation)+1)/nrofpermutations);
            p_median = ((sum(r(2)>=median_perm_correlation)+1)/nrofpermutations);
            p_mean = ((sum(r(2)>=mean_perm_correlation)+1)/nrofpermutations);
            disp(['Comparing ' num2str(c) ' - ' num2str(cc) ': r = ' num2str(r(2)) ', p_min = ' num2str(p_min) ', p_median = ' ...
                num2str(p_median) ', p_mean = ' num2str(p_mean)]);
            
        end
        
    end
    
end



