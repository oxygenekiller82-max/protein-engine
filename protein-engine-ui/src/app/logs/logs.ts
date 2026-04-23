import { Component,OnInit,inject,ChangeDetectorRef } from '@angular/core';
import { CommonModule } from '@angular/common';
import { HttpClient } from '@angular/common/http';
import { Router } from '@angular/router';
import { SelectionService } from '../services/selection.service';
//oh that's how u import it yh yh 
interface SequenceDTO {
  id: number;
  peptideChain: string;
  targetLength: number;
  isBiological: boolean;
  createdAt: string;
  masseCible: number;
  echelleKyteDoolittle: number;
  indiceAliphatique: number;
  bindingAffinity: number;
}

@Component({
  selector: 'app-logs',
  templateUrl: './logs.html',
  standalone: true,                           
  imports: [CommonModule],
  styleUrls: ['./logs.css']

  
})


 export class LogsComponent implements OnInit  {

  private http = inject(HttpClient);
  private router = inject(Router);
  private selectionService = inject(SelectionService);

  allSequences: any[] = [];
  selectedItems: any[] = [];
  private cdr = inject(ChangeDetectorRef);

 ngOnInit(){
  //my-sequenvces get
  this.http.get<SequenceDTO[]>('http://localhost:8081/api/PFA/my-sequences').subscribe({
    next: (data)=>{
      //all for the table 
      this.allSequences=data.map(s=>({ //MAP IS GENIUS HERE
        id: s.id,
        peptideChain: s.peptideChain,
        length: s.targetLength,
        isBiological: s.isBiological,
        creationTime: s.createdAt || new Date().toLocaleString(),
        status: s.peptideChain.includes('API') ? 'Failed' : 'Success',
        mass: s.masseCible,       
        hydro: s.echelleKyteDoolittle,
        stability: s.indiceAliphatique,
        binding: s.bindingAffinity
      }));
      this.cdr.detectChanges() //WAKE THE F UP ANGUALR
    }
  });
  
}

toggleSelection(seq:any){
  //????
  if (this.selectedItems.find(i => i.id === seq.id)) {
    this.selectedItems = this.selectedItems.filter(i => i.id !== seq.id);
  } else if (this.selectedItems.length < 2) {
    this.selectedItems.push(seq);
  }

  if (this.selectedItems.length === 2) {
    this.selectionService.setSelection(this.selectedItems);
    this.router.navigate(['/compare']);
  }

}
  
}