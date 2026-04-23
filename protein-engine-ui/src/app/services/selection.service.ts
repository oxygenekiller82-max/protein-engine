//select 2 rows from logs -> redirect to comapre with them there!!! 
//quick comapre! 
import { Injectable } from '@angular/core';

@Injectable({
    providedIn: 'root'
})
export class SelectionService {
    private selectedForCompare: any[] = [];//araray is fine

    setSelection(sequences: any[]) {
        this.selectedForCompare = [...sequences];
    }

    getSelection(): any[] {
        const data = [...this.selectedForCompare];
        this.selectedForCompare = []; // Clearing it... ofc 
        return data;
    }
}