import { Component, OnInit } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { Chart, registerables } from 'chart.js';

Chart.register(...registerables);

@Component({
  selector: 'app-compare',
  standalone: true,
  imports: [CommonModule, FormsModule],
  templateUrl: './compare.html',
  styleUrls: ['./compare.css']
})
export class CompareComponent implements OnInit {

  sequenceName: string = '';
  selectedSequences: string[] = [];
  storedSequences: string[] = [
    'Sequence_Alpha',
    'Sequence_Beta',
    'Sequence_Gamma',
    'Sequence_Delta'
  ];
  chartsRendered: boolean = false;

  ngOnInit(): void {}

  addSequence(seq: string): void {
    if (this.selectedSequences.length < 2 && !this.selectedSequences.includes(seq)) {
      this.selectedSequences.push(seq);
      this.sequenceName = this.selectedSequences.join(', ');
    }
  }

  removeSequence(seq: string): void {
    this.selectedSequences = this.selectedSequences.filter(s => s !== seq);
    this.sequenceName = this.selectedSequences.join(', ');
  }

  compare(): void {
    if (this.selectedSequences.length < 2) {
      alert('Veuillez sélectionner 2 séquences.');
      return;
    }
    this.chartsRendered = false;
    setTimeout(() => {
      this.chartsRendered = true;
      setTimeout(() => this.renderCharts(), 100);
    }, 0);
  }

  renderCharts(): void {
    const seq1 = this.selectedSequences[0];
    const seq2 = this.selectedSequences[1];
    const labels = ['Val 1', 'Val 2', 'Val 3', 'Val 4', 'Val 5'];
    const datasets = [
      { label: 'Hydrophobicity' },
      { label: 'Mass' },
      { label: 'Stability' },
      { label: 'Binding' },
    ];
    datasets.forEach((ds, i) => {
      const ctx = document.getElementById('chart' + i) as HTMLCanvasElement;
      if (!ctx) return;
      const existingChart = Chart.getChart(ctx);
      if (existingChart) existingChart.destroy();
      new Chart(ctx, {
        type: 'bar',
        data: {
          labels,
          datasets: [
            {
              label: seq1,
              data: labels.map(() => Math.floor(Math.random() * 100)),
              backgroundColor: 'rgba(0,0,0,0.7)',
            },
            {
              label: seq2,
              data: labels.map(() => Math.floor(Math.random() * 100)),
              backgroundColor: 'rgba(0,0,0,0.2)',
              borderColor: 'black',
              borderWidth: 1,
            }
          ]
        },
        options: {
          responsive: true,
          plugins: {
            legend: { position: 'top' },
            title: {
              display: true,
              text: ds.label + ' — ' + seq1 + ' vs ' + seq2
            }
          }
        }
      });
    });
  }
}