# Wavegraph

Wavegraph is a graph-based clustering for GW search with coherent Waveburst

## Instructions for graph production

### Template signal generation and preconditioning

* For CBC sources
```
wgsignals --input-file data/cbc_tmplbank_test.xml.gz \
          --output-file data/cbc_template_test.hdf5 \
		  --analysis-sample-rate 1024 \
		  --source-type cbc \
		  --waveform-approximant SEOBNRv2 \
		  --low-freq-cutoff 32 \
		  --wvf-gen-sample-rate 2048 \
		  --log debug
```

* Generic waveforms
```
wgsignals --input-file data/generic_tmplbank_test.hdf5 \
          --output-file data/generic_template_test.hdf5 \
		  --analysis-sample-rate 1024 \
		  --source-type generic \
		  --log debug
```

The templates can be displayed using `display_signals file.hdf5`

### Cluster computation

* For CBC sources
 * Using the matching pursuit algorithm:
```
wgclusters --input-file data/cbc_template_test.hdf5
           --output-file data/clusters_cbc_template_test.hdf5
		   --method 'matchingpursuit'
		   --segment-duration 16
		   --wvf-alignment 'right'
		   --min-scale 3
		   --max-scale 7
		   --num-time-shift 0
		   --reject-pixels-at-zero-freq True
		   --ordering 'frequency' 'time'
		   --log debug
```
 * Using the iterative hard thresholding algorithm:
 * Using Increasing Normalised Iterative Hard Thresholding:
```
wgclusters --input-file data/cbc_template_test.hdf5 \
           --output-file data/clusters_cbc_template_test.hdf5 \
		   --method 'iniht' \
		   --segment-duration 16 \
		   --wvf-alignment 'right' \
		   --min-scale 3 \
		   --max-scale 7 \
		   --num-time-shift 0 \
		   --reject-pixels-at-zero-freq False \
		   --ordering 'frequency' 'time' \
		   --log debug \
           --approx-error 0.05 \
           --cv_threshold 1e-7 \
           --c 1e-2\
           --kappa 1.5\
           --start_pixels 40 \
           --step_pixels 1 \
           --verbose False
```
* Using Constrained Basis Pursuit Denoising (L1 minimisation):
```
wgclusters --input-file data/cbc_template_test.hdf5 \
           --output-file data/clusters_cbc_template_test.hdf5 \
    	   --method 'cbpdn' \
		   --segment-duration 16 \
		   --wvf-alignment 'right' \
		   --min-scale 3 \
		   --max-scale 7 \
		   --num-time-shift 0 \
		   --reject-pixels-at-zero-freq False \
		   --ordering "frequency" "time" \
		   --log debug \
           --approx-error 0.05 \
           --cv_threshold 1e-9 \
           --cv_diff 1e-3 \
           --gamma 1e-2 \
           --kappa 10 \
           --bpdn_solver ista \
           --verbose False
```
* For generic waveforms
```
wgclusters --input-file data/generic_template_test.hdf5
           --output-file data/clusters_generic_template_test.hdf5 \
    	   --method 'matchingpursuit' \
		   --segment-duration 16 \
		   --wvf-alignment 'right' \
		   --min-scale 3 \
		   --max-scale 7 \
		   --num-time-shift 0 \
		   --reject-pixels-at-zero-freq True \
		   --ordering 'frequency' 'time' \
		   --log debug \
```

The clusters can be displayed using `display_cluster file.hdf5`
			     
## Graph for longburst waveform model

* Template signal computation:
```
wgsignals   --input-file data/example_longburst.hdf5 \
            --output-file data/example_longburst_template.hdf5 \
		    --analysis-sample-rate 4096 \
		    --source-type generic \
		    --log debug
```
* Cluster computation:
```
wgclusters  --output-file data/clusters_example_longburst.hdf5 \
		    --method 'matchingpursuit' \
		    --segment-duration 256 \
		    --wvf-alignment 'right' \
		    --min-scale 11 \
		    --max-scale 11 \
		    --num-time-shift 0 \
		    --reject-pixels-at-zero-freq True \
		    --ordering 'frequency' 'time' \
		    --log debug \
		    data/example_longburst_template.hdf5
```
* Graph generation:
```
wggraph     --input-files data/clusters_example_longburst.hdf5 \
            --output-file data/graph_example_longburst.txt \
		    --log debug
```

## Parallel computing on Condor queues

* Cluster generation:
```
wgclusters_pipe --input-files files*.hdf5 --argument-file data/demo_args.txt --sim-dir sim
```
* For CBC sources:
```
wggraph --input-files data/clusters_cbc_template_test.hdf5 \
		--output-file data/graph_cbc_template_test.txt \
		--log debug
```
* For generic waveforms:
```
wggraph --input-files data/clusters_generic_template_test.hdf5 \
        --output-file data/graph_generic_template_test.txt \
        --log debug
```

The graph can be displayed using `display_graph file.hdf5`
