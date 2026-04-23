import { Component, OnInit ,ChangeDetectorRef } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { Chart, registerables } from 'chart.js';
import { HttpClient } from '@angular/common/http';
import { inject } from '@angular/core';

import { SelectionService } from '../services/selection.service';

Chart.register(...registerables);

interface ProteinSequence {
  id: number;
  peptideChain: string;
  targetLength: number;
}

interface CompareResult {
  target_hydro: number[];
  res_hydro: number[];
  target_stability: number[];
  res_stability: number[];
  target_binding: number[];
  res_binding: number[];
  target_mass: number[];
  res_mass: number[];
  hydro_sim: number;
  stability_sim: number;
  molecular_mass_sim: number;
  binding_sim: number;
}

@Component({
  selector: 'app-compare',
  standalone: true,
  imports: [CommonModule, FormsModule],
  templateUrl: './compare.html',
  styleUrls: ['./compare.css']
})
export class CompareComponent implements OnInit {
  private http = inject(HttpClient);
  private selectionService = inject(SelectionService);
  private cdr = inject(ChangeDetectorRef);

  sequenceName: string = '';
  chartsRendered: boolean = false;

  userSavedSequences: ProteinSequence[] = [];
  selectedSequences: ProteinSequence[] = [];
  primarySequence: ProteinSequence | null = null;

  // similarity scores to display
  hydroSim: number = 0;
  stabilitySim: number = 0;
  massSim: number = 0;
  bindingSim: number = 0;

  //manually give AA
  manualSequence: string="";

  //i think will override selected sequence 2 with it if this is not ""
  //IT WORKED

  errorMessage: string ="";




  ngOnInit(): void {
    //ARE THE SEQUENCES PASSED FROM THE LOGS PAGE ?
    const preSelected = this.selectionService.getSelection();
    if (preSelected.length === 2) {
      this.selectedSequences = preSelected.map(s => ({
        id: s.id,
        peptideChain: s.peptideChain,
        targetLength: s.targetLength
      }));
      this.updateSequenceNameDisplay();
    }
      this.loadUserSequences();
  }

  loadUserSequences(): void {
    this.http.get<any[]>('http://localhost:8081/api/PFA/my-sequences')
      .subscribe({
        next: (data) => {
          this.userSavedSequences = data
          .filter(s => !s.peptideChain.includes('API')) //IF the peptide string even includes "API" it's a failed one -> no need for comparing.. because theer's
          //nothing to compare!!
            .map(s => ({
            id: s.id,
            peptideChain: s.peptideChain,
            targetLength: s.targetLength
          }));
          this.cdr.detectChanges(); //ANGULAR WAKE THE F UP
        },
        
        error: (err) =>{
           this.errorMessage="Failed to load sequences";
        }
        
        
      });
  }

  onPrimarySelect(): void {
    if (this.primarySequence) {
      this.addSequence(this.primarySequence);
      this.primarySequence = null; // reset dropdown
    }
  }

  addSequence(seq: ProteinSequence): void {
    //new checks for text area manual 
    if(this.manualSequence.trim() !=="" && this.selectedSequences.length>0){
      this.errorMessage = "clear the manual sequence input to select a second sequence.";
      return;
    }
    if (this.selectedSequences.length < 2 && !this.selectedSequences.some(s => s.id === seq.id)) {
      this.selectedSequences.push(seq);
      this.updateSequenceNameDisplay();

    }else if(this.selectedSequences.length>=2){
        this.errorMessage = "a maximum of 2 sequences can be added";
    }
    
  }

  removeSequence(seq: ProteinSequence): void {
    this.selectedSequences = this.selectedSequences.filter(s => s.id !== seq.id);
    this.chartsRendered = false; // hide charts if selection changes
    this.updateSequenceNameDisplay();
    this.errorMessage = "";
  }

  private updateSequenceNameDisplay(): void {
    this.sequenceName = this.selectedSequences.map(s => `Seq #${s.id}`).join(' vs ');
  }

  compare(): void {
    if (this.selectedSequences.length===0 ) {
      this.errorMessage ="Please select at least one sequence from the dropdown menu";
      return;
    }

    if (this.selectedSequences.length===0 && this.manualSequence.trim()==="" ) {
      this.errorMessage ="Comparison required 2 sequences.";
      return;
    }

    //drop down = 2 + manual 
    if (this.selectedSequences.length === 2 && this.manualSequence.trim() !== "") {
      this.errorMessage = "Clear the manual input or remove a sequence — you already have 2 selected.";
      return;
    }

    //just 1 selected + no manula

    if (this.selectedSequences.length === 1 && this.manualSequence.trim() === "") {
      this.errorMessage = "Second sequence is missing — select one or type it manually.";
      return;
    }


    let seq_two = this.manualSequence!== "" ?
     this.manualSequence :
     (this.selectedSequences[1].peptideChain||"");

    
     if (this.selectedSequences.length<2 && this.manualSequence.trim()==="") {
      this.errorMessage = "Second sequence is missing (select one or type it).";
      return;

    }
    //LENGTH must 
    const payload = {
      target_seq: this.selectedSequences[0].peptideChain,
      generated_seq: seq_two
      
    };

    this.http.post<CompareResult>('http://localhost:8081/api/PFA/compare', payload)
      .subscribe({
        next: (data: CompareResult) => {
          // store similarity scores
          this.hydroSim = data.hydro_sim;
          this.stabilitySim = data.stability_sim;
          this.massSim = data.molecular_mass_sim;
          this.bindingSim = data.binding_sim;

          this.chartsRendered = true;
          setTimeout(() => this.renderCharts(data), 100);
        },
        error: (err) =>{
          this.errorMessage="Server Error: " + (err.error?.message || "Server connection failed.");
         console.error('Comparison failed', err)
        }
      });

  }

  renderCharts(data: CompareResult): void {
    const seq1Label = `Seq #${this.selectedSequences[0].id}`;
    const seq2Label = this.manualSequence.trim() !== "" 
  ? "Manual Input" 
  : `Seq #${this.selectedSequences[1].id}`;

    const chartConfigs = [
      { canvasId: 'hydroChart',     label: 'Hydropathy',       seq1Values: data.target_hydro,     seq2Values: data.res_hydro },
      { canvasId: 'stabilityChart', label: 'Stability',        seq1Values: data.target_stability, seq2Values: data.res_stability },
      { canvasId: 'bindingChart',   label: 'Binding Affinity', seq1Values: data.target_binding,   seq2Values: data.res_binding },
      { canvasId: 'massChart',      label: 'Molecular Mass',   seq1Values: data.target_mass,      seq2Values: data.res_mass },
    ];

    chartConfigs.forEach(config => {
      const canvas = document.getElementById(config.canvasId) as HTMLCanvasElement;
      if (!canvas) return;

      // destroy existing chart to avoid overlap on re-compare
      const existing = Chart.getChart(canvas);
      if (existing) existing.destroy();

      const labels = config.seq1Values.map((_: number, i: number) => `${i + 1}`);

      new Chart(canvas, {
        type: 'line',
        data: {
          labels,
          datasets: [
            {
              label: `${seq1Label} (target)`,
              data: config.seq1Values,
              borderColor: '#4f8ef7',
              backgroundColor: 'rgba(79,142,247,0.1)',
              borderWidth: 2,
              pointRadius: 0,
              tension: 0.3,
            },
            {
              label: `${seq2Label} (generated)`,
              data: config.seq2Values,
              borderColor: '#f97316',
              backgroundColor: 'rgba(249,115,22,0.1)',
              borderWidth: 2,
              borderDash: [5, 5],
              pointRadius: 0,
              tension: 0.3,
            },
          ],
        },
        options: {
          responsive: true,
          plugins: {
            legend: { position: 'top' },
            title: {
              display: true,
              text: `${config.label} Comparison`,
            },
          },
          scales: {
            x: { title: { display: true, text: 'Amino Acid Position' } },
            y: { title: { display: true, text: config.label } },
          },
        },
      });
    });
  }
}